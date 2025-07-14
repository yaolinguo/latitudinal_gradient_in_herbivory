# Data_zr_herbivory-latitude.csv

This dataset contains effect sizes (`zr`) quantifying the strength and direction of latitudinal herbivory gradients across native and non-native plant species.

## Column descriptions

| Column name              | Description |
|--------------------------|-------------|
| `Reference`              | Abbreviated citation for the original study |
| `Journal_name`           | Name of the journal where the study was published |
| `Publication_year`       | Year the study was published |
| `Experiment_year`        | Year when the original herbivory experiment or field data were collected |
| `IF`                     | Journal impact factor (5-year IF from Web of Science) |
| `Data_details`           | Notes about data extraction (e.g., whether derived from figures, tables, raw data) |
| `Variable_name`          | Description of the herbivory type|
| `Variable_identity`      | Data were based on herbivory intensity (e.g., percentage of leaf area removed) or on herbivore abundance, and were used as proxies for herbivory pressure in our study |
| `Plant_species`          | Scientific name of the plant species involved |
| `Plant_status`           | Whether the plant is native or non-native in the study location (`native` / `non-native`) |
| `Herbivore_feeding_guild`| Type of herbivore interacting with the plant (`All folivores`, `defoliator`, `Grazers`, `Gallers`, `Miners`, `Seed feeders`, `Sap feeders`, `Stem feeders`) |
| `cor`                    | Raw Pearson correlation coefficient between herbivory and latitude |
| `n`                      | Number of observations (e.g., sites or populations) used to calculate the correlation |
| `zr`                     | Fisher’s z-transformed correlation coefficient |
| `var`                    | Variance of the `zr` effect size |

## Notes

- Negative `zr` values indicate increasing herbivory at lower latitudes.
- Used in REML models to compare overall gradients and plant status-specific patterns.




# Data_zr_herbivory-bioclimate.csv

This dataset contains Fisher’s z‑transformed effect sizes (zr), calculated from correlations between herbivory intensity and 19 bioclimatic variables (Bio1–Bio19). It was used to evaluate how herbivory patterns respond to climatic gradients across latitude.

## Column descriptions

- `Data_zr_herbivory-bioclimate.csv` contains the same columns as `Data_zr_herbivory-latitude.csv`, with one additional column:
  
| Column name       | Description |
|-------------------|-------------|
| `Bio_var`         | Indicates the specific bioclimatic variable used to calculate each effect size. Bioclimatic predictors include annual mean temperature (Bio 1), mean diurnal range (Bio 2), isothermality (Bio 3), temperature seasonality (Bio 4), maximum temperature of the warmest month (Bio 5), minimum temperature of the coldest month (Bio 6), annual temperature range (Bio 7), mean temperature of the wettest quarter (Bio 8), mean temperature of the driest quarter (Bio 9), mean temperature of the warmest quarter (Bio 10), mean temperature of the coldest quarter (Bio 11), annual precipitation (Bio 12), precipitation of the wettest month (Bio 13), precipitation of the driest month (Bio 14), precipitation seasonality (Bio 15), precipitation of the wettest quarter (Bio 16), precipitation of the driest quarter (Bio 17), precipitation of the warmest quarter (Bio 18), and precipitation of the coldest quarter (Bio 19). |

## Notes

- Bioclimatic variables follow WorldClim naming (e.g., `Bio1` = Annual Mean Temperature, `Bio12` = Annual Precipitation).
- This dataset supports moderator analysis of climatic drivers on herbivory.
