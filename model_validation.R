library(dismo)
library(gbm)
library(pROC)
library(dplyr)

file_path_train <- "train_data.csv"
train_df <- read.csv(file_path_train)

target <- names(train_df)[1]
all_vars <- names(train_df)[2:8]

k_folds <- 5
set.seed(123)
folds <- sample(rep(1:k_folds, length.out = nrow(train_df)))

run_cv <- function(data, predictors, target, folds,
                   tree.complexity = 5, learning.rate = 0.005, bag.fraction = 0.5) {
  auc_vals <- numeric(max(folds))
  
  x_idx <- match(predictors, names(data))
  y_idx <- match(target, names(data))
  
  for (k in 1:max(folds)) {
    train_data <- data[folds != k, , drop = FALSE]
    test_data  <- data[folds == k, , drop = FALSE]
    
    model <- gbm.step(
      data = train_data,
      gbm.x = x_idx,
      gbm.y = y_idx,
      family = "bernoulli",
      tree.complexity = tree.complexity,
      learning.rate = learning.rate,
      bag.fraction = bag.fraction,
      verbose = FALSE
    )
    
    preds <- predict.gbm(
      model,
      test_data,
      n.trees = model$gbm.call$best.trees,
      type = "response"
    )
    
    roc_obj <- roc(test_data[[target]], preds)
    auc_vals[k] <- as.numeric(auc(roc_obj))
  }
  
  auc_vals
}

results <- data.frame()

full_auc <- run_cv(train_df, all_vars, target, folds)
full_mean <- mean(full_auc)
full_sd <- sd(full_auc)

results <- rbind(results, data.frame(
  Model = "Full model (all variables)",
  Removed = "-",
  Fold1 = round(full_auc[1], 2),
  Fold2 = round(full_auc[2], 2),
  Fold3 = round(full_auc[3], 2),
  Fold4 = round(full_auc[4], 2),
  Fold5 = round(full_auc[5], 2),
  Mean = round(full_mean, 2),
  SD = round(full_sd, 2),
  Delta_AUC = 0.00
))

for (removed in all_vars) {
  sel_vars <- setdiff(all_vars, removed)
  auc_vals <- run_cv(train_df, sel_vars, target, folds)
  
  results <- rbind(results, data.frame(
    Model = "Leave-one-out",
    Removed = removed,
    Fold1 = round(auc_vals[1], 2),
    Fold2 = round(auc_vals[2], 2),
    Fold3 = round(auc_vals[3], 2),
    Fold4 = round(auc_vals[4], 2),
    Fold5 = round(auc_vals[5], 2),
    Mean = round(mean(auc_vals), 2),
    SD = round(sd(auc_vals), 2),
    Delta_AUC = round(mean(auc_vals) - full_mean, 2)
  ))
  
  cat("Finished removing:", removed, "- Mean AUC =", round(mean(auc_vals), 3), "\n")
}

results <- results %>%
  mutate(`Mean (±SD)` = sprintf("%.2f (±%.02f)", Mean, SD)) %>%
  select(Model, Removed, Fold1:Fold5, `Mean (±SD)`, Delta_AUC)

print(results)

write.csv(
  results,
  "crossvalidation_performance.csv",
  row.names = FALSE
)
