# Light Exposure and Alzheimer’s Disease (Mediation Analysis)

## Overview
This repository contains R scripts for mediation analyses examining the association between light exposure and risk of Alzheimer’s disease (AD), using UK Biobank data.

The primary aim is to investigate whether biological and behavioral factors mediate the association between light exposure and AD risk.

## Study Objectives
To evaluate whether the effect of light exposure on Alzheimer’s disease is mediated by:

- Brain structural measures
- Circadian rest–activity rhythm (CRAR) parameters
- Serum vitamin D levels

## Data Source
UK Biobank:
- Accelerometer-derived light exposure metrics
- Brain MRI imaging-derived phenotypes (IDPs)
- CRAR parameters
- Blood biomarkers (e.g., vitamin D)
- Dementia outcomes and covariates

Note: Raw UK Biobank data are not included in this repository.

## Repository Structure

- `Mediation_Brain.R` – Mediation analysis using brain structural measures
- `Mediation_CRAR.R` – Mediation analysis using circadian rhythm parameters
- `Mediation_VitaminD.R` – Mediation analysis using serum vitamin D

## Statistical Methods

- Cox proportional hazards models
- Mediation analysis (e.g., R package `mediation` or equivalent framework)
- Covariate adjustment (age, sex, education, lifestyle, comorbidities, etc.)
- Multiple testing correction where applicable

## Software Requirements

- R (≥ 4.2)
- Required packages may include:
  - survival
  - mediation
  - dplyr
  - broom
  - car

## Reproducibility
This repository contains analysis scripts only.  
Researchers with approved UK Biobank access can reproduce results using the same analytical pipeline.

## Author
Wei Wang
