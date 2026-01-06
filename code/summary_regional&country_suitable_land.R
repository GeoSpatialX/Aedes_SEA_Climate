suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

base_dir  <- "output"
gcm_list  <- c("ACCESS-CM2","EC-Earth3","MIROC6","MPI-ESM1-2-HR","NorESM2-MM","UKESM1-0-LL")
ssps      <- c("SSP126","SSP245","SSP585")
years_all <- c(2000, 2050, 2080)

use_iqr <- FALSE      # TRUE: IQR band (25–75%); FALSE: min–max
top_n   <- 20         # number of countries to show in the country-change plot

baseline_ae_path <- "ae_baseline.tif"
baseline_al_path <- "al_baseline.tif"

thr_ae <- ae_threshold
thr_al <- al_threshold

countries_shp <- "sea.shp"
countries <- st_read(countries_shp, quiet = TRUE) |> st_make_valid()

out_dir <- file.path(base_dir, "summary_outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

find_pred <- function(gcm, species, ssp, year) {
  d <- file.path(base_dir, gcm, paste0(species, "output"), "new_edition")
  if (!dir.exists(d)) return(character(0))
  pat <- paste0("^", species, ".*_", ssp, "_", year, "\\.tif$")
  list.files(d, pattern = pat, full.names = TRUE)
}

share_over_thr <- function(r, thr, land_mask = NULL) {
  if (!is.null(land_mask)) r <- mask(r, land_mask)
  hit <- r >= thr
  ok  <- !is.na(r)
  num <- global(hit, "sum", na.rm = TRUE)[1, 1]
  den <- global(ok,  "sum", na.rm = TRUE)[1, 1]
  100 * as.numeric(num / den)
}

pick_country_name <- function(sf_obj) {
  nm <- names(sf_obj)
  cand <- c("name","NAME","admin","ADMIN","NAME_EN","NAME_0","name_0","COUNTRY")
  use <- cand[cand %in% nm][1]
  if (is.na(use)) stop("Country name field not found in the shapefile.")
  sf_obj[[use]]
}

read_future <- function(gcm, species, ssp, year, ref_r, mask_ref) {
  f <- find_pred(gcm, species, ssp, year)
  if (!length(f)) return(NULL)
  r <- rast(f[1])
  r <- project(r, ref_r, method = "bilinear")
  r <- mask(r, mask_ref)
  r
}

share_by_country <- function(r, thr, countries_sf, land_mask = NULL) {
  if (!is.null(land_mask)) r <- mask(r, land_mask)
  hit <- r >= thr; names(hit) <- "hit"
  ok  <- !is.na(r); names(ok) <- "ok"
  ext <- terra::extract(c(hit, ok), vect(st_transform(countries_sf, crs(r))), df = TRUE)
  
  ext |>
    filter(!is.na(ok)) |>
    group_by(ID) |>
    summarise(
      share = 100 * sum(hit, na.rm = TRUE) / sum(ok, na.rm = TRUE),
      .groups = "drop"
    )
}

summ_band <- function(x, use_iqr = FALSE) {
  tibble(
    mean = mean(x, na.rm = TRUE),
    med  = median(x, na.rm = TRUE),
    lo   = if (use_iqr) as.numeric(quantile(x, 0.25, na.rm = TRUE)) else min(x, na.rm = TRUE),
    hi   = if (use_iqr) as.numeric(quantile(x, 0.75, na.rm = TRUE)) else max(x, na.rm = TRUE),
    n    = length(x)
  )
}

summarise_delta <- function(df, tag, use_iqr = FALSE) {
  df |>
    group_by(country) |>
    summarise(
      mu  = mean(delta, na.rm = TRUE),
      lo  = if (use_iqr) quantile(delta, 0.25, na.rm = TRUE) else min(delta, na.rm = TRUE),
      hi  = if (use_iqr) quantile(delta, 0.75, na.rm = TRUE) else max(delta, na.rm = TRUE),
      n   = n(),
      .groups = "drop"
    ) |>
    mutate(kind = tag)
}

run_regional_summary <- function(species, baseline_path, thr, use_iqr = FALSE) {
  message("Regional summary: ", toupper(species))
  
  r0 <- rast(baseline_path)
  share2000 <- share_over_thr(r0, thr)
  
  rows <- list(
    tibble(year = 2000, ssp = "Baseline", mean = share2000, med = share2000, lo = share2000, hi = share2000, n = 1)
  )
  
  for (yr in years_all[years_all != 2000]) {
    for (ssp in ssps) {
      vals <- c()
      for (gcm in gcm_list) {
        f <- find_pred(gcm, species, ssp, yr)
        if (!length(f)) next
        rv <- rast(f[1])
        vals <- c(vals, share_over_thr(rv, thr))
      }
      if (length(vals)) {
        rows[[length(rows) + 1]] <- tibble(
          year = yr, ssp = ssp,
          mean = mean(vals), med = median(vals),
          lo   = if (use_iqr) quantile(vals, 0.25) else min(vals),
          hi   = if (use_iqr) quantile(vals, 0.75) else max(vals),
          n    = length(vals)
        )
      }
    }
  }
  
  whole_df <- bind_rows(rows)
  
  pal <- c("Baseline"="grey60","SSP126"="#76B7B2","SSP245"="#F28E2B","SSP585"="#E15759")
  pd  <- position_dodge(width = 0.5)
  
  p <- ggplot(whole_df, aes(x = factor(year), y = mean, color = ssp)) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, position = pd, linewidth = 0.5) +
    geom_point(size = 3.2, position = pd) +
    scale_color_manual(values = pal, name = NULL) +
    labs(x = NULL, y = "Suitable land area (%)") +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "top",
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.text.y  = element_text(size = 18),
      axis.text.x  = element_text(size = 18),
      legend.text  = element_text(size = 16)
    )
  
  ggsave(file.path(out_dir, paste0(species, "_regional_suitable_land_by_ssp.png")),
         p, width = 7.5, height = 4.6, dpi = 300)
  
  list(data = whole_df, plot = p)
}

country_change_two_deltas <- function(species, baseline_path, thr,
                                      countries_sf, use_iqr = FALSE,
                                      top_n = 20, fixed_levels = NULL) {
  message("Country-level change: ", toupper(species))
  
  r0 <- rast(baseline_path)
  mask_ref <- !is.na(r0)
  
  countries_sf$cn_name <- pick_country_name(countries_sf)
  
  share0_by_cty <- share_by_country(r0, thr, countries_sf, land_mask = mask_ref) |>
    mutate(country = countries_sf$cn_name[ID]) |>
    select(country, share0 = share)
  
  shares <- list()
  for (yr in c(2050, 2080)) {
    for (ssp in ssps) {
      for (gcm in gcm_list) {
        rv <- read_future(gcm, species, ssp, yr, r0, mask_ref)
        if (is.null(rv)) next
        s_cty <- share_by_country(rv, thr, countries_sf, land_mask = mask_ref) |>
          mutate(country = countries_sf$cn_name[ID], year = yr, ssp = ssp, gcm = gcm) |>
          select(country, year, ssp, gcm, share)
        shares[[length(shares) + 1]] <- s_cty
      }
    }
  }
  shares <- bind_rows(shares)
  if (!nrow(shares)) stop("No valid prediction files were found.")
  
  delta_50_00 <- shares |>
    filter(year == 2050) |>
    left_join(share0_by_cty, by = "country") |>
    mutate(delta = share - share0) |>
    select(country, ssp, gcm, delta)
  
  wide_50_80 <- shares |>
    filter(year %in% c(2050, 2080)) |>
    pivot_wider(names_from = year, values_from = share) |>
    filter(!is.na(`2050`), !is.na(`2080`)) |>
    mutate(delta = `2080` - `2050`) |>
    select(country, ssp, gcm, delta)
  
  sum_50_00 <- summarise_delta(delta_50_00, "2050–2000", use_iqr = use_iqr)
  sum_80_50 <- summarise_delta(wide_50_80,   "2080–2050", use_iqr = use_iqr)
  both <- bind_rows(sum_50_00, sum_80_50)
  
  if (is.null(fixed_levels)) {
    ord <- sum_50_00 |>
      arrange(desc(mu)) |>
      slice_head(n = top_n)
    keep_levels <- ord$country
  } else {
    keep_levels <- fixed_levels
  }
  
  plot_df <- both |>
    filter(country %in% keep_levels) |>
    mutate(
      country = factor(country, levels = rev(keep_levels)),
      kind    = factor(kind, levels = c("2080–2050", "2050–2000"))
    )
  
  pd <- position_dodge(width = 0.6)
  p <- ggplot(plot_df, aes(x = country, y = mu, color = kind)) +
    geom_errorbar(aes(ymin = lo, ymax = hi), position = pd, width = 0.5, linewidth = 1) +
    geom_point(position = pd, size = 3) +
    coord_flip() +
    scale_color_manual(
      breaks = c("2080–2050", "2050–2000"),
      values = c("2080–2050" = "#0d95ce", "2050–2000" = "#fbcb1f"),
      name = NULL
    ) +
    labs(x = NULL, y = "Change in suitable land area (%)") +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "top",
      legend.text  = element_text(size = 16),
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.text.y  = element_text(size = 18),
      axis.text.x  = element_text(size = 18),
      axis.title.x = element_text(size = 18)
    )
  
  ggsave(file.path(out_dir, paste0(species, "_country_change_two_deltas.png")),
         p, width = 8.2, height = 6.6, dpi = 600)
  
  write.csv(both, file.path(out_dir, paste0(species, "_country_change_two_deltas_table.csv")),
            row.names = FALSE)
  
  list(plot = p, levels = keep_levels, table = both)
}

res_reg_ae <- run_regional_summary("ae", baseline_ae_path, thr_ae, use_iqr = use_iqr)
res_reg_al <- run_regional_summary("al", baseline_al_path, thr_al, use_iqr = use_iqr)

res_cty_ae <- country_change_two_deltas(
  "ae", baseline_ae_path, thr_ae, countries,
  use_iqr = use_iqr, top_n = top_n
)

res_cty_al <- country_change_two_deltas(
  "al", baseline_al_path, thr_al, countries,
  use_iqr = use_iqr, top_n = top_n,
  fixed_levels = res_cty_ae$levels
)
