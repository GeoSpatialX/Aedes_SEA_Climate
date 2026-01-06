# Load required libraries
library(dismo)
library(pROC)
library(gbm)
library(terra)

# --------------------------
# 1) Train GBM model (AE)
# --------------------------

file_path_ae_train <- "ae_train_data"
AE_train <- read.csv(file_path_ae_train)

ae_model <- gbm.step(
  data = AE_train,
  gbm.x = 2:8,
  gbm.y = 1,
  family = "bernoulli",
  tree.complexity = 5,
  learning.rate = 0.005,
  bag.fraction = 0.5
)

# --------------------------
# 2) Evaluate on test data (AE)
# --------------------------

file_path_ae_test <- "ae_test_data"
AE_test <- read.csv(file_path_ae_test)

preds_ae <- predict(
  ae_model,
  AE_test,
  n.trees = ae_model$gbm.call$best.trees,
  type = "response"
)

roc_ae <- roc(AE_test$PRESENCE, preds_ae)
best_threshold_ae <- coords(
  roc_ae,
  x = "best",
  ret = "threshold",
  best.method = "closest.topleft"
)$threshold

plot(
  roc_ae,
  col = "orange",
  main = paste("AE AUC-ROC =", round(auc(roc_ae), 2)),
  ylim = c(0, 1),
  xlim = c(1, 0),
  xaxs = "i",
  yaxs = "i"
)
abline(a = 0, b = 1, lty = 2, col = "red")

# --------------------------
# 3) Predict suitability on rasters (AE)
# --------------------------
# This block can be used for baseline and future projections.
# Replace the input TIFF with the corresponding baseline/future raster stack
# (same band order and same predictor names).

ae_tif <- "ae_Combined_7bands.tif"
AE_grids <- rast(ae_tif)

band_names <- c(
  "Temperature_Suitability",
  "Precipitation_of_Driest_Month",
  "Precipitation_of_Wettest_Month",
  "Relative_Humidity",
  "Urban_Fraction",
  "DEM",
  "POP"
)
names(AE_grids) <- band_names

p_ae <- predict(
  AE_grids,
  ae_model,
  type = "response",
  fun = predict.gbm,
  n.trees = ae_model$gbm.call$best.trees,
  na.rm = TRUE
)

output_continuous_ae <- "ae_environmental_suitability.tif"
writeRaster(p_ae, output_continuous_ae, overwrite = TRUE)

# --------------------------
# 4) AL follows the same workflow
# --------------------------
# Train/evaluate AL using the AL training/testing CSVs, then replace the raster stack with
# the AL baseline/future predictor TIFFs (same band names/order).