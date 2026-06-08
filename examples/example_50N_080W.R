# =========================================================
# EXAMPLE: SINGLE TILE RUN (50N_080W)
# =========================================================

library(terra)
library(caret)
library(ranger)
library(gstat)
library(ggplot2)
library(viridisLite)

# -------------------------
# LOAD FUNCTIONS
# -------------------------

source("R/io_utils.R")
source("R/raster_utils.R")
source("R/modeling_utils.R")
source("R/mapping_utils.R")
source("scripts/run_tile.R")

# -------------------------
# RUN SINGLE TILE
# -------------------------

result <- run_tile(
  tile_id = "50N_080W",
  tile_file = "./Data/GLAD_Tiles/50N_080W.tif",
  reference_file = "./Data/dAGB_3km.tif",
  output_dir = "./Outputs/50N_080W"
)

# -------------------------
# QUICK VISUAL CHECK
# -------------------------

plot_raster(result$final_map, "Final Biomass Map (50N_080W)")

cat("\nExample completed successfully\n")
