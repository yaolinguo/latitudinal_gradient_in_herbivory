# Data_zr_herbivory-latitude.csv

This UTF‑8–encoded CSV holds Fisher **z‑transformed effect sizes (`zr`)** quantifying the relationship between herbivory and latitude for both native and non‑native plant species.

## Column descriptions

| Column name              | Description |
|--------------------------|-------------|
| `Reference`              | Abbreviated citation for the original study (First author and year). |
| `Journal_name`           | Name of the journal where the study was published. |
| `Publication_year`       | Year of publication. |
| `Experiment_year`        | Year(s) in which herbivory data were collected.|
| `IF`                     | Five‑year journal impact factor (Web of Science, 2024 edition). |
| `Data_details`           | Notes about data extraction (e.g., whether derived from figures, tables, raw data) |
| `Variable_name`          | Short description of the herbivory type|
| `Variable_identity`      | Categorical flag: `intensity` (direct damage such as % leaf‑area removed) or `abundance` (herbivore counts used as a proxy for herbivory pressure). |
| `Plant_species`          | Scientific name of the plant species involved |
| `Plant_status`           | Whether the plant is native or non-native in the study location (`native` / `non-native`) |
| `Herbivore_feeding_guild`| Type of herbivore feeding guild: `All folivores`, `defoliators`, `Grazers`, `Gallers`, `Miners`, `Seed feeders`, `Sap feeders`, `Stem feeders`. |
| `cor`                    | Raw Pearson correlation coefficient between herbivory and latitude |
| `n`                      | Number of observations (e.g., sites or populations) used to calculate the correlation |
| `zr`                     | Fisher’s z-transformed correlation coefficient |
| `var`                    | Variance of the `zr` effect size |

## Notes

- Negative `zr` values indicate stronger herbivory toward lower latitudes; positive values indicate the opposite.
- These effect sizes feed into REML meta‑analytic models comparing overall and status‑specific latitudinal gradients.

---

# Data_zr_herbivory-bioclimate.csv

A companion UTF‑8 CSV that expands the latitude dataset with **WorldClim 2.1 bioclimatic predictors**.  Each row is the correlation between herbivory and one bioclimatic variable for a given plant–study combination.

## Column descriptions

- **All columns are identical to _Data_zr_herbivory-latitude.csv_**, plus one additional identifier:
  
| Column name       | Description |
|-------------------|-------------|
| `Bio_var`         | Code for the bioclimatic variable used in the correlation:<br>• `Bio1` – Annual mean temperature<br>• `Bio2` – Mean diurnal range<br>• `Bio3` – Isothermality<br>• `Bio4` – Temperature seasonality<br>• `Bio5` – Max temperature of warmest month<br>• `Bio6` – Min temperature of coldest month<br>• `Bio7` – Annual temperature range<br>• `Bio8` – Mean temp. of wettest quarter<br>• `Bio9` – Mean temp. of driest quarter<br>• `Bio10` – Mean temp. of warmest quarter<br>• `Bio11` – Mean temp. of coldest quarter<br>• `Bio12` – Annual precipitation<br>• `Bio13` – Precipitation of wettest month<br>• `Bio14` – Precipitation of driest month<br>• `Bio15` – Precipitation seasonality<br>• `Bio16` – Precipitation of wettest quarter<br>• `Bio17` – Precipitation of driest quarter<br>• `Bio18` – Precipitation of warmest quarter<br>• `Bio19` – Precipitation of coldest quarter |

## Notes

- Units follow WorldClim conventions.  
- Negative `zr` implies herbivory increases as the climate variable decreases (and vice‑versa).  
- Used as moderators in mixed‑effects REML models examining bioclimatic effects on herbivory across latitude.

