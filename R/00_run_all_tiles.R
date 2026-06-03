# =========================
# 00_run_all_tiles.R
# GLOBAL 228 TILE PIPELINE
# =========================

library(terra)
library(dplyr)

tiles <- read.csv("tiles_list.csv")

source("R/01_data_prep.R")
source("R/02_rf_model.R")
source("R/03_rf_prediction.R")
source("R/04_kriging.R")
source("R/05_final_mapping.R")

# -------------------------
# OUTPUT BASE DIR
# -------------------------
output_base <- "Results"

if (!dir.exists(output_base)) {
  dir.create(output_base)
}

# -------------------------
# LOOP OVER TILES
# -------------------------

for (i in 1:nrow(tiles)) {
  
  fname <- tiles$tile_id[i]
  
  cat("\n=====================================\n")
  cat("Processing tile:", fname, "\n")
  cat("=====================================\n")
  
  try({
    
    # Set global variable for scripts
    assign("fname", fname, envir = .GlobalEnv)
    
    # -------------------------
    # RUN PIPELINE
    # -------------------------
    
    source("R/01_data_prep.R")
    source("R/02_rf_model.R")
    source("R/03_rf_prediction.R")
    source("R/04_kriging.R")
    source("R/05_final_mapping.R")
    
    cat("SUCCESS:", fname, "\n")
    
  }, silent = FALSE)
  
  cat("DONE TILE:", fname, "\n\n")
}