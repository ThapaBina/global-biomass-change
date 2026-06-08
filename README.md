# Global Biomass Change Mapping (RF + Kriging Framework)

## Overview

This repository provides a reproducible workflow for high-resolution global biomass mapping using a hybrid machine learning and geostatistical framework.

The method combines:

- Random Forest modeling
- Residual kriging (geostatistical correction)
- Multi-resolution raster fusion
- Mass-preservation scaling

---

## Key Features

- Tile-based processing (200+ global tiles)
- Fully reproducible pipeline
- Parallel execution support
- RF + kriging hybrid modeling
- 30 m and 3 km multi-resolution outputs
- Publication-ready figures and diagnostics

---

## Repository Structure

---

## Example Usage

### Run single tile

```r
source("scripts/run_tile.R")

run_tile(
  tile_id = "50N_080W",
  tile_file = "Data/GLAD_Tiles/50N_080W.tif",
  reference_file = "Data/dAGB_3km.tif",
  output_dir = "Outputs/50N_080W"
)
