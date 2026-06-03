# =========================
# 02_rf_model.R
# =========================

library(caret)
library(ranger)

train_df <- readRDS("outputs/train_df.rds")

set.seed(42)

ctrl <- trainControl(
  method = "cv",
  number = 5,
  verboseIter = TRUE,
  savePredictions = "final"
)

grid <- expand.grid(
  mtry = 2,
  splitrule = "variance",
  min.node.size = 5
)

rf_model <- train(
  dAGB ~ tch_mean + tch_sd + tch_q10 + tch_q90,
  data = train_df,
  method = "ranger",
  trControl = ctrl,
  tuneGrid = grid,
  num.trees = 1000,
  importance = "permutation",
  quantreg = TRUE
)

saveRDS(rf_model, "outputs/rf_model.rds")