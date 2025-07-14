# Data_zr_herbivory-latitude.csv

This dataset contains effect sizes (`zr`) quantifying the strength and direction of latitudinal herbivory gradients across native and non-native plant species.

## Purpose

Used to evaluate whether latitudinal herbivory patterns differ in magnitude and direction between native and non-native plants.

## Column descriptions

| Column name       | Description |
|-------------------|-------------|
| `Study_ID`        | Unique identifier for each reference/study (usually corresponds to the reference or citation) |
| `Plant_species`   | Scientific name of the focal plant species observed in the study |
| `Plant_status`    | Indicates whether the plant is native or non-native in the sampled region (`native` / `non-native`) |
| `Feeding_guild`   | Herbivore functional group (e.g., `defoliator`, `sap-feeder`, `galler`, etc.) |
| `zr`              | Fisher’s z-transformed Pearson correlation coefficient between herbivory and latitude |
| `Var_zr`          | Sampling variance of the `zr` value |
| `N_sites`         | Number of geographic sampling sites used in the correlation |
| `Latitude_min`    | Minimum latitude of the sampling range |
| `Latitude_max`    | Maximum latitude of the sampling range |
| `Year`            | Year the data were collected or published |
| `Reference`       | Full or abbreviated citation of the source study |

## Notes

- Negative `zr` values indicate increasing herbivory at lower latitudes.
- Used in REML models to compare overall trends and status-specific patterns.
