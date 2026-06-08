# Global Biomass Change Mapping Workflow

## Overview

This repository implements a hybrid machine learning + geostatistical framework for mapping aboveground biomass (AGB) using:

- Random Forest (RF)
- Residual kriging (ATPK-style correction)
- Mass-preservation scaling
- Multi-resolution raster fusion

---

## Methodological Pipeline

### Step 1: Input Data

- 30 m GLAD-derived raster predictors
- 3 km reference biomass (dAGB)

---

### Step 2: Feature Aggregation

Fine resolution (30 m) → aggregated to 3 km:

- Mean
- Standard deviation
- Quantiles (10%, 90%)

---

### Step 3: Random Forest Modeling

Model:
dAGB ~ tch_mean + tch_sd + tch_q10 + tch_q90

Outputs:
- RF prediction (3 km)
- Variable importance
- Cross-validation metrics

---

### Step 4: Residual Calculation
Residual = Observed - RF Prediction

---

### Step 5: Variogram Modeling

- Empirical variogram computed on residuals
- Spherical model fitted

---

### Step 6: Kriging

- Ordinary kriging applied to RF residuals
- Produces spatial correction surface

---

### Step 7: RF + Kriging Fusion
Final = RF + Kriging Residuals

---

### Step 8: Downscaling to 30 m

- 3 km correction surface resampled to 30 m
- Combined with RF prediction

---

### Step 9: Mass Preservation Correction

Ensures consistency with reference biomass:
Correction = Reference / Aggregated Prediction

Applied multiplicatively at 30 m resolution.

---

## Outputs

Each tile produces:

- RF prediction (3 km)
- Kriging residuals (3 km)
- Final corrected biomass (30 m)
- Variance maps
- Figures for validation

---

## Scaling

The framework is designed for:

- ~220 tiles
- parallel execution
- cluster computing

---

## Software Dependencies

- R (>= 4.2)
- terra
- caret
- ranger
- gstat
- ggplot2
- future.apply











