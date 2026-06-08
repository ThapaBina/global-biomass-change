# =========================================================
# RUN ALL TILES (GLOBAL BIOMASS CHANGE PIPELINE)
# =========================================================

suppressPackageStartupMessages({
  library(parallel)
  library(future.apply)
})

source("R/io_utils.R")
source("scripts/run_tile.R")

# -------------------------
# CONFIGURATION
# -------------------------

input_dir <- "./Data/GLAD_Tiles"
ref_file  <- "./Data/dAGB_3km.tif"
output_dir <- "./Outputs"

ensure_dir(output_dir)

tile_files <- list.files(input_dir, pattern = "\\.tif$", full.names = TRUE)

tile_ids <- tools::file_path_sans_ext(basename(tile_files))

# -------------------------
# PARALLEL SETUP
# -------------------------

plan <- future::plan(future::multisession, workers = max(1, parallel::detectCores() - 1))

# -------------------------
# SAFE WRAPPER
# -------------------------

run_safe <- function(tile_file, tile_id) {

  out_tile_dir <- file.path(output_dir, tile_id)
  ensure_dir(out_tile_dir)

  cat("\n====================================\n")
  cat("STARTING TILE:", tile_id, "\n")
  cat("====================================\n")

  tryCatch({

    run_tile(
      tile_id = tile_id,
      tile_file = tile_file,
      reference_file = ref_file,
      output_dir = out_tile_dir
    )

    cat("SUCCESS:", tile_id, "\n")

  }, error = function(e) {

    cat("FAILED:", tile_id, "\n")
    cat("ERROR:", conditionMessage(e), "\n")

    writeLines(
      paste0(tile_id, " | ERROR: ", conditionMessage(e)),
      file.path(output_dir, "error_log.txt"),
      append = TRUE
    )
  })
}

# -------------------------
# RUN ALL TILES
# -------------------------

future.apply::future_mapply(
  FUN = run_safe,
  tile_file = tile_files,
  tile_id = tile_ids,
  SIMPLIFY = FALSE
)

cat("\nALL TILES COMPLETED\n")
