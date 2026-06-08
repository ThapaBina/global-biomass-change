# =========================
# 01_data_prep.R
# =========================

library(terra)
library(dplyr)

source("R/utils.R")

setwd("D:/Example")

fname <- "50N_080W"

infile <- file.path("Data/GLAD_Tiles", paste0(fname, ".tif"))

ba_3km <- rast("Data/dAGB_3km.tif")
ba_30m <- rast(infile)

# Crop
ba_3km <- crop(ba_3km, ext(ba_30m))
names(ba_3km) <- "dAGB"

# aggregation factor
fact <- round(res(ba_3km)[1] / res(ba_30m)[1])

# predictors
tch_mean <- aggregate(ba_30m, fact, mean, na.rm = TRUE)
tch_sd   <- aggregate(ba_30m, fact, sd, na.rm = TRUE)
tch_q    <- aggregate(ba_30m, fact, q_fun)

names(tch_q) <- c("q10", "q90")

tch_mean <- resample(tch_mean, ba_3km)
tch_sd   <- resample(tch_sd, ba_3km)
tch_q    <- resample(tch_q, ba_3km)

pred_3km <- c(tch_mean, tch_sd, tch_q)
names(pred_3km) <- c("tch_mean","tch_sd","tch_q10","tch_q90")

train_df <- as.data.frame(c(ba_3km, pred_3km), xy = TRUE, na.rm = TRUE)

saveRDS(train_df, "outputs/train_df.rds")
saveRDS(pred_3km, "outputs/pred_3km.rds")
saveRDS(ba_3km, "outputs/ba_3km.rds")