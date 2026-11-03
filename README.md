# Marine Fish Species Distribution Modeling (SDM)

## Acanthochromis polyacanthus

This repository presents an R-based species distribution modeling workflow investigating the potential distribution of *Acanthochromis polyacanthus* across three climatic periods:

- **Past:** Last Glacial Maximum (LGM; ~21,000 years ago)
- **Present:** Bio-ORACLE v3 environmental conditions
- **Future:** SSP5-8.5 climate scenario (2090–2100)

### Methods

Species occurrence records were obtained using GBIF. Environmental predictors include sea surface temperature, salinity, and bathymetry, using Bio-ORACLE v3 and paleo-MARSPEC datasets.

Species distribution models were developed in R using a binomial generalized linear model (GLM).

### Repository Contents

- `scripts/` – R scripts for GBIF occurrence extraction and SDM analysis
- `figures/` – Maps of modeled species distributions
- `results/` – Summary statistics, environmental predictor correlations, and comparisons of suitable habitat area

### Research Objectives

To investigate potential shifts in suitable marine habitats across past, present, and future climatic conditions, with particular emphasis on the Philippine marine environment.

### Data Sources

- GBIF: https://www.gbif.org/
- Bio-ORACLE: https://www.bio-oracle.org/
- MARSPEC: https://lifewatch.github.io/sdmpredictors/

### Notes

Model predictions represent potential environmental suitability rather than confirmed species occurrence. Historical and future projections should be interpreted in light of environmental data resolution and modeling assumptions.