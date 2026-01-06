# Load packages
library(dismo)
library(pROC)
library(gbm)
library(raster)
library(ggplot2)
library(viridis)
library(dplyr)

# 1) Fit the AE model
file_path <- "ae_train.csv"
Ae_train <- read.csv(file_path)

head(Ae_train)
sum(is.na(Ae_train))

ae_model <- gbm.step(
  data = Ae_train,
  gbm.x = 2:8,
  gbm.y = 1,
  family = "bernoulli",
  tree.complexity = 5,
  learning.rate = 0.01,
  bag.fraction = 0.5
)

# 2) Sensitivity analysis setup
input_directory  <- "ae_avg_future_suitability"
output_directory <- "sensitivity_output"

if (!dir.exists(output_directory)) dir.create(output_directory, recursive = TRUE)

tif_files <- list.files(
  path = input_directory,
  pattern = "\\.tif$",
  full.names = TRUE
)

sensitivity_results <- data.frame()
sensitivity_rasters <- list()

replacement_tif <- stack("ae_combined_7bands_baseline.tif")


# 3) Run sensitivity analysis
for (tif_file_path in tif_files) {
  
  ssp  <- gsub(".*_(SSP[0-9]+)_(\\d{4})\\.tif$", "\\1", basename(tif_file_path))
  year <- gsub(".*_(\\d{4})\\.tif$", "\\1", basename(tif_file_path))
  
  ae_grids <- brick(tif_file_path)
  names(ae_grids) <- c(
    "Temperature_Suitability",
    "Precipitation_of_Driest_Month",
    "Precipitation_of_Wettest_Month",
    "Relative_Humidity",
    "Urban_Fraction",
    "DEM",
    "POP"
  )
  
  ae_grids <- projectRaster(ae_grids, crs = crs(replacement_tif))
  ae_grids <- resample(ae_grids, replacement_tif, method = "bilinear")
  ae_grids <- crop(ae_grids, replacement_tif)
  ae_grids <- mask(ae_grids, replacement_tif)
  
  baseline_prediction <- predict(
    ae_grids, ae_model,
    n.trees = ae_model$gbm.call$best.trees,
    type = "response"
  )
  
  for (i in 1:7) {
    variable <- names(ae_grids)[i]
    
    ae_grids_fixed <- ae_grids
    ae_grids_fixed[[variable]] <- replacement_tif[[i]]
    
    fixed_prediction <- predict(
      ae_grids_fixed, ae_model,
      n.trees = ae_model$gbm.call$best.trees,
      type = "response"
    )
    
    sensitivity <- (fixed_prediction - baseline_prediction) / baseline_prediction * 100
    
    sensitivity_results <- rbind(
      sensitivity_results,
      data.frame(
        Variable = variable,
        SSP = ssp,
        Year = year,
        Mean_Sensitivity = cellStats(sensitivity, stat = "mean", na.rm = TRUE)
      )
    )
    
    sensitivity_rasters[[paste0(variable, "_", ssp, "_", year)]] <- sensitivity
    
    cat("Processed:", variable, ssp, year, "\n")
  }
}

# 4) Plot by SSP and year
for (year in unique(sensitivity_results$Year)) {
  for (ssp in unique(sensitivity_results$SSP)) {
    
    year_data <- subset(sensitivity_results, Year == year & SSP == ssp)
    if (nrow(year_data) == 0) next
    
    year_data <- year_data |>
      group_by(Variable) |>
      summarise(Mean_Sensitivity = mean(Mean_Sensitivity)) |>
      ungroup()
    
    year_data$Plot_Value <- year_data$Mean_Sensitivity
    year_data$Plot_Value[year_data$Variable == "Temperature_Suitability"] <-
      year_data$Plot_Value[year_data$Variable == "Temperature_Suitability"] / 30
    
    p <- ggplot(year_data, aes(x = Plot_Value, y = Variable, fill = Variable)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7) +
      geom_bar(stat = "identity", width = 0.5) +
      geom_text(
        aes(
          x = Plot_Value + ifelse(Plot_Value > 0, 0.15, -0.15),
          label = round(Plot_Value, 2)
        ),
        hjust = ifelse(year_data$Plot_Value > 0, 0, 1),
        size = 7
      ) +
      scale_fill_manual(values = c("#041F4A", "#08306B", "#08519C", "#2171B5",
                                   "#6BAED6", "#C6DBEF", "#E8F1FA")) +
      coord_cartesian(xlim = c(-6, 1)) +
      theme_classic() +
      theme(
        legend.position = "none",
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 18)
      ) +
      labs(x = "Sensitivity", y = "Variable")
    
    ggsave(
      filename = file.path(output_directory, paste0("ae_sensitivity_", ssp, "_", year, ".png")),
      plot = p, width = 12, height = 7, dpi = 600
    )
  }
}

# The same workflow applies to AL by switching the training CSV, the input raster folder,
# and the baseline replacement stack to the AL equivalents.
