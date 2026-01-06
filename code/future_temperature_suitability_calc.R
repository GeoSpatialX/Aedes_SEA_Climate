# This script builds on the temperature suitability framework described in:
# Brady, O. J. et al. (2014). Global temperature constraints on Aedes aegypti
# and Ae. albopictus persistence and competence for dengue virus transmission.
# Parasites & Vectors, 7, 338.
#
# The original code has been adapted to generate temperature suitability
# surfaces for Aedes aegypti and Aedes albopictus, with modifications for
# Southeast Asia and for application to baseline and future climate scenarios.

rm(list = ls())

library(snowfall)
library(Rcpp)
library(tempsuitcalc)
library(raster)
library(tidyverse)

# Set species: "aegypti" or "albopictus"
Species <- "aegypti"

# Load weekly minimum and maximum temperature data
load("weekly_max_temperature.RData")
load("weekly_min_temperature.RData")

if (Species == "aegypti") {
  load("tempmatAe2.RData")
  Tempmat <- tempmatAe2
  weekmin <- tmin_weekly_matrix_all
  weekmax <- tmax_weekly_matrix_all
} else if (Species == "albopictus") {
  load("tempmatAl2.RData")
  Tempmat <- tempmatAl2
  weekmin <- tmin_weekly_matrix_all
  weekmax <- tmax_weekly_matrix_all
} else {
  stop("Species not recognised. Use 'aegypti' or 'albopictus'.")
}

# Extend weekly temperature series to account for burn-in and burn-out weeks
weekmin <- cbind(weekmin[, 48:52], weekmin, weekmin[, 4:8])
weekmax <- cbind(weekmax[, 48:52], weekmax, weekmax[, 4:8])

# Identify grid cells with valid temperature information
row_sum <- rowSums(weekmin)
valid_cells <- row_sum >= 1
cell_index <- seq_len(nrow(weekmin))

x <- data.frame(cell_id = cell_index[valid_cells])
x <- cbind(x, weekmin[valid_cells, ], weekmax[valid_cells, ])

# Directory containing processed temperature suitability CSV files
tempsuit_dir <- "C:/Users/lls/Desktop/global_new/ae"
tempsuit_files <- list.files(
  tempsuit_dir,
  pattern = "_processed.csv$",
  full.names = TRUE
)

# Template raster defining spatial resolution and extent
template_raster <- raster("weekly_temp.tif")

# Convert temperature suitability matrices to raster stacks
for (tempsuit_file in tempsuit_files) {
  
  tempsuit <- read_csv(tempsuit_file)
  tempsuit <- as.matrix(tempsuit[, 3:108])
  
  if (nrow(tempsuit) != nrow(x)) {
    stop("Mismatch between temperature suitability matrix and valid cell index.")
  }
  
  final_matrix <- matrix(
    NA,
    nrow = ncell(template_raster),
    ncol = ncol(tempsuit)
  )
  
  valid_rows <- !is.na(x$cell_id)
  cell_ids <- x$cell_id[valid_rows]
  tempsuit_valid <- tempsuit[valid_rows, ]
  
  final_matrix[cell_ids, ] <- tempsuit_valid
  
  raster_stack <- stack()
  for (i in seq_len(ncol(final_matrix))) {
    r <- raster(template_raster)
    r[] <- final_matrix[, i]
    raster_stack <- addLayer(raster_stack, r)
  }
  
  model_name <- basename(dirname(dirname(tempsuit_file)))
  ssp <- str_extract(basename(tempsuit_file), "ssp[0-9]+")
  
  output_file <- paste0(
    tempsuit_dir, "/",
    model_name, "_",
    ssp, "_time_", Species, "_tempsuit.tif"
  )
  
  writeRaster(
    raster_stack,
    filename = output_file,
    format = "GTiff",
    overwrite = TRUE
  )
  
  cat("Processed and saved:", output_file, "\n")
}
