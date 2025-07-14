#########################################################################################################
### Guo et al., 2025. Herbivory increases towards lower latitudes in native but not introduced plants ###
#########################################################################################################

library(knitr)
library(tidyverse)
library(brms)
library(ape)
library(coda)
library(modelr)
library(gridExtra)
library(pBrackets)
library(RColorBrewer)
library(performance)
library(phytools)
library(kableExtra)
library(tidybayes)
library(formattable)
library(grid)
library(dplyr)
library(ggplot2)
library(base)
library(bayesplot)
library(metafor)
library(ggbeeswarm)
library(pander)
library(raster)
library(maps)
library(mapdata)
library(rworldmap)
library(readxl)
library(utils)
library(reshape2)
library(ggpubr)
library(stringr)
library(ggstance)
library(ggridges) 
library(sp)
library(raster)
library(geodata)
library(multcompView)
library(multcomp)
library(rotl)
library(ape)
library(psych)
library(meta)
library(devtools)
library(processx)
library(metafor)
library(gghalves)
library(readxl)
library(ape)
library(phytools)
library(writexl)
library(dplyr)
library(tidybayes)
library(tidybayes)
library(ggplot2)
library(dplyr)
library(ggdist)
library(forcats)
library(scales)
library(geodata)
library(terra)
library(sf)
library(dplyr)
library(purrr)
library(metafor)
library(purrr)
library(broom.mixed)

### Data cleaning ###
setwd("/Users/yaolin/Desktop/My papers/Manuscripts/2025 - Guo - BioRxiv - Meta for herbivory pattern/Submission to PNAS - 0527")
data <- read_excel("Raw_data_zr - 0705.xlsx")
selected_cols <- c("Reference", "Journal_name", "Publication_year", "Experiment_year", "IF",
                   "Data_details", "Variable_name", "Variable_identity", "Plant_species",
                   "Study_region", "Study_region_detail", "Plant_status", "Native_region",
                   "Native_region_detail", "Hemisphere_EW", "Hemisphere_SN", "Clime",
                   "Methods", "Invasive_year", "Invasive_duration", "Herbivore_feeding_guild",
                   "Herbivore_species", "Herbivore_detail", "Herbivore_diet_breadth",
                   "Herbivore_diet_breadth_detail", "Longitude", "Latitude", "Herbivory")
data <- data[, selected_cols]

outdir <- file.path(getwd(), "worldclim_data")
bio <- worldclim_global(var = "bio", res = 10, path = outdir) 
data <- data %>% mutate(across(c(Longitude, Latitude), as.numeric)) %>% 
  filter(!is.na(Longitude) & !is.na(Latitude))

pts_sf <- st_as_sf(data, coords = c("Longitude", "Latitude"), crs = 4326) 
pts_vect <- vect(pts_sf)
bio_vals <- extract(bio, pts_vect)
data_bioclim <- bind_cols(data, extract(bio, pts_vect, ID = FALSE) %>% as.data.frame())
write.csv(data_bioclim, "data_bioclim.csv", row.names = FALSE)

data_bioclim <- data_bioclim %>% rename_with(~ sub("^wc2\\.1_10m_bio_", "Bio_", .x), starts_with("wc2.1_10m_bio_"))
group_vars <- c("Reference", "Journal_name", "Publication_year", "Experiment_year", "IF",
                "Data_details", "Variable_name", "Variable_identity", "Plant_species",
                "Plant_status", "Herbivore_feeding_guild")

# Bio_1
data_bioclim_zr_Bio_1 <- data_bioclim %>%
  filter(!is.na(Bio_1) & !is.na(Herbivory)) %>%
  group_by(across(all_of(group_vars))) %>%
  summarise(cor = cor(Bio_1, Herbivory, use = "complete.obs", method = "pearson"),
            n   = sum(!is.na(Bio_1) & !is.na(Herbivory)), .groups = "drop") %>%
  mutate(zr = atanh(cor), var = ifelse(n == 3, 1, 1/(n-3)))

rm1_Bio_1 <- rma.mv(yi = zr, V = var, mods = ~ 1,
                    random = ~ 1 | Plant_species + Reference,
                    method = "REML", data = data_bioclim_zr_Bio_1)

rm2_Bio_1 <- rma.mv(yi = zr, V = var, mods = ~ Plant_status - 1,
                    random = ~ 1 | Plant_species + Reference,
                    method = "REML", data = data_bioclim_zr_Bio_1)

# Bio_2
data_bioclim_zr_Bio_2 <- data_bioclim %>%
  filter(!is.na(Bio_2) & !is.na(Herbivory)) %>%
  group_by(across(all_of(group_vars))) %>%
  summarise(cor = cor(Bio_2, Herbivory, use = "complete.obs", method = "pearson"),
            n   = sum(!is.na(Bio_2) & !is.na(Herbivory)), .groups = "drop") %>%
  mutate(zr = atanh(cor), var = ifelse(n == 3, 1, 1/(n-3)))

rm1_Bio_2 <- rma.mv(yi = zr, V = var, mods = ~ 1,
                    random = ~ 1 | Plant_species + Reference,
                    method = "REML", data = data_bioclim_zr_Bio_2)

rm2_Bio_2 <- rma.mv(yi = zr, V = var, mods = ~ Plant_status - 1,
                    random = ~ 1 | Plant_species + Reference,
                    method = "REML", data = data_bioclim_zr_Bio_2)

# Bio_3
data_bioclim_zr_Bio_3 <- data_bioclim %>%
  filter(!is.na(Bio_3) & !is.na(Herbivory)) %>%
  group_by(across(all_of(group_vars))) %>%
  summarise(cor = cor(Bio_3, Herbivory, use = "complete.obs", method = "pearson"),
            n   = sum(!is.na(Bio_3) & !is.na(Herbivory)), .groups = "drop") %>%
  mutate(zr = atanh(cor), var = ifelse(n == 3, 1, 1/(n-3)))

rm1_Bio_3 <- rma.mv(yi = zr, V = var, mods = ~ 1,
                    random = ~ 1 | Plant_species + Reference,
                    method = "REML", data = data_bioclim_zr_Bio_3)

rm2_Bio_3 <- rma.mv(yi = zr, V = var, mods = ~ Plant_status - 1,
                    random = ~ 1 | Plant_species + Reference,
                    method = "REML", data = data_bioclim_zr_Bio_3)

# Bio_4
data_bioclim_zr_Bio_4 <- data_bioclim %>%
  filter(!is.na(Bio_4) & !is.na(Herbivory)) %>%
  group_by(across(all_of(group_vars))) %>%
  summarise(cor = cor(Bio_4, Herbivory, use = "complete.obs", method = "pearson"),
            n   = sum(!is.na(Bio_4) & !is.na(Herbivory)), .groups = "drop") %>%
  mutate(zr = atanh(cor), var = ifelse(n == 3, 1, 1/(n-3)))

rm1_Bio_4 <- rma.mv(yi = zr, V = var, mods = ~ 1,
                    random = ~ 1 | Plant_species + Reference,
                    method = "REML", data = data_bioclim_zr_Bio_4)

rm2_Bio_4 <- rma.mv(yi = zr, V = var, mods = ~ Plant_status - 1,
                    random = ~ 1 | Plant_species + Reference,
                    method = "REML", data = data_bioclim_zr_Bio_4)

# Bio_5
data_bioclim_zr_Bio_5 <- data_bioclim %>%
  filter(!is.na(Bio_5) & !is.na(Herbivory)) %>%
  group_by(across(all_of(group_vars))) %>%
  summarise(cor = cor(Bio_5, Herbivory, use = "complete.obs", method = "pearson"),
            n   = sum(!is.na(Bio_5) & !is.na(Herbivory)), .groups = "drop") %>%
  mutate(zr = atanh(cor), var = ifelse(n == 3, 1, 1/(n-3)))

rm1_Bio_5 <- rma.mv(yi = zr, V = var, mods = ~ 1,
                    random = ~ 1 | Plant_species + Reference,
                    method = "REML", data = data_bioclim_zr_Bio_5)

rm2_Bio_5 <- rma.mv(yi = zr, V = var, mods = ~ Plant_status - 1,
                    random = ~ 1 | Plant_species + Reference,
                    method = "REML", data = data_bioclim_zr_Bio_5)

# Bio_6
data_bioclim_zr_Bio_6 <- data_bioclim %>%
  filter(!is.na(Bio_6) & !is.na(Herbivory)) %>%
  group_by(across(all_of(group_vars))) %>%
  summarise(cor = cor(Bio_6, Herbivory, use = "complete.obs", method = "pearson"),
            n   = sum(!is.na(Bio_6) & !is.na(Herbivory)), .groups = "drop") %>%
  mutate(zr = atanh(cor), var = ifelse(n == 3, 1, 1/(n-3)))

rm1_Bio_6 <- rma.mv(yi = zr, V = var, mods = ~ 1,
                    random = ~ 1 | Plant_species + Reference,
                    method = "REML", data = data_bioclim_zr_Bio_6)

rm2_Bio_6 <- rma.mv(yi = zr, V = var, mods = ~ Plant_status - 1,
                    random = ~ 1 | Plant_species + Reference,
                    method = "REML", data = data_bioclim_zr_Bio_6)

# Bio_7
data_bioclim_zr_Bio_7 <- data_bioclim %>%
  filter(!is.na(Bio_7) & !is.na(Herbivory)) %>%
  group_by(across(all_of(group_vars))) %>%
  summarise(cor = cor(Bio_7, Herbivory, use = "complete.obs", method = "pearson"),
            n   = sum(!is.na(Bio_7) & !is.na(Herbivory)), .groups = "drop") %>%
  mutate(zr = atanh(cor), var = ifelse(n == 3, 1, 1/(n-3)))

rm1_Bio_7 <- rma.mv(yi = zr, V = var, mods = ~ 1,
                    random = ~ 1 | Plant_species + Reference,
                    method = "REML", data = data_bioclim_zr_Bio_7)

rm2_Bio_7 <- rma.mv(yi = zr, V = var, mods = ~ Plant_status - 1,
                    random = ~ 1 | Plant_species + Reference,
                    method = "REML", data = data_bioclim_zr_Bio_7)

# Bio_8
data_bioclim_zr_Bio_8 <- data_bioclim %>%
  filter(!is.na(Bio_8) & !is.na(Herbivory)) %>%
  group_by(across(all_of(group_vars))) %>%
  summarise(cor = cor(Bio_8, Herbivory, use = "complete.obs", method = "pearson"),
            n   = sum(!is.na(Bio_8) & !is.na(Herbivory)), .groups = "drop") %>%
  mutate(zr = atanh(cor), var = ifelse(n == 3, 1, 1/(n-3)))

rm1_Bio_8 <- rma.mv(yi = zr, V = var, mods = ~ 1,
                    random = ~ 1 | Plant_species + Reference,
                    method = "REML", data = data_bioclim_zr_Bio_8)

rm2_Bio_8 <- rma.mv(yi = zr, V = var, mods = ~ Plant_status - 1,
                    random = ~ 1 | Plant_species + Reference,
                    method = "REML", data = data_bioclim_zr_Bio_8)

# Bio_9
data_bioclim_zr_Bio_9 <- data_bioclim %>%
  filter(!is.na(Bio_9) & !is.na(Herbivory)) %>%
  group_by(across(all_of(group_vars))) %>%
  summarise(cor = cor(Bio_9, Herbivory, use = "complete.obs", method = "pearson"),
            n   = sum(!is.na(Bio_9) & !is.na(Herbivory)), .groups = "drop") %>%
  mutate(zr = atanh(cor), var = ifelse(n == 3, 1, 1/(n-3)))

rm1_Bio_9 <- rma.mv(yi = zr, V = var, mods = ~ 1,
                    random = ~ 1 | Plant_species + Reference,
                    method = "REML", data = data_bioclim_zr_Bio_9)

rm2_Bio_9 <- rma.mv(yi = zr, V = var, mods = ~ Plant_status - 1,
                    random = ~ 1 | Plant_species + Reference,
                    method = "REML", data = data_bioclim_zr_Bio_9)

# Bio_10
data_bioclim_zr_Bio_10 <- data_bioclim %>%
  filter(!is.na(Bio_10) & !is.na(Herbivory)) %>%
  group_by(across(all_of(group_vars))) %>%
  summarise(cor = cor(Bio_10, Herbivory, use = "complete.obs", method = "pearson"),
            n   = sum(!is.na(Bio_10) & !is.na(Herbivory)), .groups = "drop") %>%
  mutate(zr = atanh(cor), var = ifelse(n == 3, 1, 1/(n-3)))

rm1_Bio_10 <- rma.mv(yi = zr, V = var, mods = ~ 1,
                     random = ~ 1 | Plant_species + Reference,
                     method = "REML", data = data_bioclim_zr_Bio_10)

rm2_Bio_10 <- rma.mv(yi = zr, V = var, mods = ~ Plant_status - 1,
                     random = ~ 1 | Plant_species + Reference,
                     method = "REML", data = data_bioclim_zr_Bio_10)

# Bio_11
data_bioclim_zr_Bio_11 <- data_bioclim %>%
  filter(!is.na(Bio_11) & !is.na(Herbivory)) %>%
  group_by(across(all_of(group_vars))) %>%
  summarise(cor = cor(Bio_11, Herbivory, use = "complete.obs", method = "pearson"),
            n   = sum(!is.na(Bio_11) & !is.na(Herbivory)), .groups = "drop") %>%
  mutate(zr = atanh(cor), var = ifelse(n == 3, 1, 1/(n-3)))

rm1_Bio_11 <- rma.mv(yi = zr, V = var, mods = ~ 1,
                     random = ~ 1 | Plant_species + Reference,
                     method = "REML", data = data_bioclim_zr_Bio_11)

rm2_Bio_11 <- rma.mv(yi = zr, V = var, mods = ~ Plant_status - 1,
                     random = ~ 1 | Plant_species + Reference,
                     method = "REML", data = data_bioclim_zr_Bio_11)

# Bio_12
data_bioclim_zr_Bio_12 <- data_bioclim %>%
  filter(!is.na(Bio_12) & !is.na(Herbivory)) %>%
  group_by(across(all_of(group_vars))) %>%
  summarise(cor = cor(Bio_12, Herbivory, use = "complete.obs", method = "pearson"),
            n   = sum(!is.na(Bio_12) & !is.na(Herbivory)), .groups = "drop") %>%
  mutate(zr = atanh(cor), var = ifelse(n == 3, 1, 1/(n-3)))

rm1_Bio_12 <- rma.mv(yi = zr, V = var, mods = ~ 1,
                     random = ~ 1 | Plant_species + Reference,
                     method = "REML", data = data_bioclim_zr_Bio_12)

rm2_Bio_12 <- rma.mv(yi = zr, V = var, mods = ~ Plant_status - 1,
                     random = ~ 1 | Plant_species + Reference,
                     method = "REML", data = data_bioclim_zr_Bio_12)

# Bio_13
data_bioclim_zr_Bio_13 <- data_bioclim %>%
  filter(!is.na(Bio_13) & !is.na(Herbivory)) %>%
  group_by(across(all_of(group_vars))) %>%
  summarise(cor = cor(Bio_13, Herbivory, use = "complete.obs", method = "pearson"),
            n   = sum(!is.na(Bio_13) & !is.na(Herbivory)), .groups = "drop") %>%
  mutate(zr = atanh(cor), var = ifelse(n == 3, 1, 1/(n-3)))

rm1_Bio_13 <- rma.mv(yi = zr, V = var, mods = ~ 1,
                     random = ~ 1 | Plant_species + Reference,
                     method = "REML", data = data_bioclim_zr_Bio_13)

rm2_Bio_13 <- rma.mv(yi = zr, V = var, mods = ~ Plant_status - 1,
                     random = ~ 1 | Plant_species + Reference,
                     method = "REML", data = data_bioclim_zr_Bio_13)

# Bio_14
data_bioclim_zr_Bio_14 <- data_bioclim %>%
  filter(!is.na(Bio_14) & !is.na(Herbivory)) %>%
  group_by(across(all_of(group_vars))) %>%
  summarise(cor = cor(Bio_14, Herbivory, use = "complete.obs", method = "pearson"),
            n   = sum(!is.na(Bio_14) & !is.na(Herbivory)), .groups = "drop") %>%
  mutate(zr = atanh(cor), var = ifelse(n == 3, 1, 1/(n-3)))

rm1_Bio_14 <- rma.mv(yi = zr, V = var, mods = ~ 1,
                     random = ~ 1 | Plant_species + Reference,
                     method = "REML", data = data_bioclim_zr_Bio_14)

rm2_Bio_14 <- rma.mv(yi = zr, V = var, mods = ~ Plant_status - 1,
                     random = ~ 1 | Plant_species + Reference,
                     method = "REML", data = data_bioclim_zr_Bio_14)

# Bio_15
data_bioclim_zr_Bio_15 <- data_bioclim %>%
  filter(!is.na(Bio_15) & !is.na(Herbivory)) %>%
  group_by(across(all_of(group_vars))) %>%
  summarise(cor = cor(Bio_15, Herbivory, use = "complete.obs", method = "pearson"),
            n   = sum(!is.na(Bio_15) & !is.na(Herbivory)), .groups = "drop") %>%
  mutate(zr = atanh(cor), var = ifelse(n == 3, 1, 1/(n-3)))

rm1_Bio_15 <- rma.mv(yi = zr, V = var, mods = ~ 1,
                     random = ~ 1 | Plant_species + Reference,
                     method = "REML", data = data_bioclim_zr_Bio_15)

rm2_Bio_15 <- rma.mv(yi = zr, V = var, mods = ~ Plant_status - 1,
                     random = ~ 1 | Plant_species + Reference,
                     method = "REML", data = data_bioclim_zr_Bio_15)

# Bio_16
data_bioclim_zr_Bio_16 <- data_bioclim %>%
  filter(!is.na(Bio_16) & !is.na(Herbivory)) %>%
  group_by(across(all_of(group_vars))) %>%
  summarise(cor = cor(Bio_16, Herbivory, use = "complete.obs", method = "pearson"),
            n   = sum(!is.na(Bio_16) & !is.na(Herbivory)), .groups = "drop") %>%
  mutate(zr = atanh(cor), var = ifelse(n == 3, 1, 1/(n-3)))

rm1_Bio_16 <- rma.mv(yi = zr, V = var, mods = ~ 1,
                     random = ~ 1 | Plant_species + Reference,
                     method = "REML", data = data_bioclim_zr_Bio_16)

rm2_Bio_16 <- rma.mv(yi = zr, V = var, mods = ~ Plant_status - 1,
                     random = ~ 1 | Plant_species + Reference,
                     method = "REML", data = data_bioclim_zr_Bio_16)

# Bio_17
data_bioclim_zr_Bio_17 <- data_bioclim %>%
  filter(!is.na(Bio_17) & !is.na(Herbivory)) %>%
  group_by(across(all_of(group_vars))) %>%
  summarise(cor = cor(Bio_17, Herbivory, use = "complete.obs", method = "pearson"),
            n   = sum(!is.na(Bio_17) & !is.na(Herbivory)), .groups = "drop") %>%
  mutate(zr = atanh(cor), var = ifelse(n == 3, 1, 1/(n-3)))

rm1_Bio_17 <- rma.mv(yi = zr, V = var, mods = ~ 1,
                     random = ~ 1 | Plant_species + Reference,
                     method = "REML", data = data_bioclim_zr_Bio_17)

rm2_Bio_17 <- rma.mv(yi = zr, V = var, mods = ~ Plant_status - 1,
                     random = ~ 1 | Plant_species + Reference,
                     method = "REML", data = data_bioclim_zr_Bio_17)

# Bio_18
data_bioclim_zr_Bio_18 <- data_bioclim %>%
  filter(!is.na(Bio_18) & !is.na(Herbivory)) %>%
  group_by(across(all_of(group_vars))) %>%
  summarise(cor = cor(Bio_18, Herbivory, use = "complete.obs", method = "pearson"),
            n   = sum(!is.na(Bio_18) & !is.na(Herbivory)), .groups = "drop") %>%
  mutate(zr = atanh(cor), var = ifelse(n == 3, 1, 1/(n-3)))

rm1_Bio_18 <- rma.mv(yi = zr, V = var, mods = ~ 1,
                     random = ~ 1 | Plant_species + Reference,
                     method = "REML", data = data_bioclim_zr_Bio_18)

rm2_Bio_18 <- rma.mv(yi = zr, V = var, mods = ~ Plant_status - 1,
                     random = ~ 1 | Plant_species + Reference,
                     method = "REML", data = data_bioclim_zr_Bio_18)

# Bio_19
data_bioclim_zr_Bio_19 <- data_bioclim %>%
  filter(!is.na(Bio_19) & !is.na(Herbivory)) %>%
  group_by(across(all_of(group_vars))) %>%
  summarise(cor = cor(Bio_19, Herbivory, use = "complete.obs", method = "pearson"),
            n   = sum(!is.na(Bio_19) & !is.na(Herbivory)), .groups = "drop") %>%
  mutate(zr = atanh(cor), var = ifelse(n == 3, 1, 1/(n-3)))

rm1_Bio_19 <- rma.mv(yi = zr, V = var, mods = ~ 1,
                     random = ~ 1 | Plant_species + Reference,
                     method = "REML", data = data_bioclim_zr_Bio_19)

rm2_Bio_19 <- rma.mv(yi = zr, V = var, mods = ~ Plant_status - 1,
                     random = ~ 1 | Plant_species + Reference,
                     method = "REML", data = data_bioclim_zr_Bio_19)

# creat dataset for bioclimatic effect sizes
bio_vars <- paste0("Bio_", 1:19)        
combined_zr <- map_dfr(bio_vars, function(bio) {
  obj_name <- paste0("data_bioclim_zr_", bio) 
  if (!exists(obj_name, envir = .GlobalEnv)) {
    message("Can not find：", obj_name, "Skip")
    return(NULL)
  }
  get(obj_name, envir = .GlobalEnv) %>% 
    mutate(Bio_var = bio)
})

write_csv(combined_zr, "data_bioclim_zr_Bio1_19_combined.csv")






# REML modeling
sample.sizes <- as.data.frame(table(data_bioclim_zr_Bio_1$Plant_status))
sample.sizes
sample.sizes <- as.data.frame(table(data_bioclim_zr_Bio_2$Plant_status))
sample.sizes
sample.sizes <- as.data.frame(table(data_bioclim_zr_Bio_3$Plant_status))
sample.sizes
sample.sizes <- as.data.frame(table(data_bioclim_zr_Bio_4$Plant_status))
sample.sizes
sample.sizes <- as.data.frame(table(data_bioclim_zr_Bio_5$Plant_status))
sample.sizes
sample.sizes <- as.data.frame(table(data_bioclim_zr_Bio_6$Plant_status))
sample.sizes
sample.sizes <- as.data.frame(table(data_bioclim_zr_Bio_7$Plant_status))
sample.sizes
sample.sizes <- as.data.frame(table(data_bioclim_zr_Bio_8$Plant_status))
sample.sizes
sample.sizes <- as.data.frame(table(data_bioclim_zr_Bio_9$Plant_status))
sample.sizes
sample.sizes <- as.data.frame(table(data_bioclim_zr_Bio_10$Plant_status))
sample.sizes
sample.sizes <- as.data.frame(table(data_bioclim_zr_Bio_11$Plant_status))
sample.sizes
sample.sizes <- as.data.frame(table(data_bioclim_zr_Bio_12$Plant_status))
sample.sizes
sample.sizes <- as.data.frame(table(data_bioclim_zr_Bio_13$Plant_status))
sample.sizes
sample.sizes <- as.data.frame(table(data_bioclim_zr_Bio_14$Plant_status))
sample.sizes
sample.sizes <- as.data.frame(table(data_bioclim_zr_Bio_15$Plant_status))
sample.sizes
sample.sizes <- as.data.frame(table(data_bioclim_zr_Bio_16$Plant_status))
sample.sizes
sample.sizes <- as.data.frame(table(data_bioclim_zr_Bio_17$Plant_status))
sample.sizes
sample.sizes <- as.data.frame(table(data_bioclim_zr_Bio_18$Plant_status))
sample.sizes
sample.sizes <- as.data.frame(table(data_bioclim_zr_Bio_19$Plant_status))
sample.sizes

rm1_Bio_1
rm2_Bio_1

rm1_Bio_2
rm2_Bio_2

rm1_Bio_3
rm2_Bio_3

rm1_Bio_4
rm2_Bio_4

rm1_Bio_5
rm2_Bio_5

rm1_Bio_6
rm2_Bio_6

rm1_Bio_7
rm2_Bio_7

rm1_Bio_8
rm2_Bio_8

rm1_Bio_9
rm2_Bio_9

rm1_Bio_10
rm2_Bio_10

rm1_Bio_11
rm2_Bio_11

rm1_Bio_12
rm2_Bio_12

rm1_Bio_13
rm2_Bio_13

rm1_Bio_14
rm2_Bio_14

rm1_Bio_15
rm2_Bio_15

rm1_Bio_16
rm2_Bio_16

rm1_Bio_17
rm2_Bio_17

rm1_Bio_18
rm2_Bio_18

rm1_Bio_19
rm2_Bio_19

model_names <- paste0(rep(c("rm1_Bio_", "rm2_Bio_"), each = 19), 1:19)
model_results <- map_dfr(model_names, function(m) {
  model <- get(m)
  tidy_result <- broom::tidy(model, conf.int = TRUE)
  tidy_result$model_name <- m
  tidy_result
})

model_results <- model_results %>%
  dplyr::select(model_name, everything())
write.csv(model_results, "model_results.csv", row.names = FALSE)
head(model_results)

# Figure 4a
bio_vars <- paste0("Bio_", 1:19)

df3 <- lapply(bio_vars, function(bio) {
  mdl <- get(paste0("rm1_", bio), envir = .GlobalEnv) 
  est <- coef(mdl)[1]                             
  se  <- sqrt(vcov(mdl)[1, 1])
  data.frame(Guild = bio,
             est   = est,
             lo    = est - 1.96 * se,
             hi    = est + 1.96 * se,
             stringsAsFactors = FALSE)
}) |> bind_rows() |>
  mutate(Guild = factor(Guild, levels = rev(bio_vars))) 

p4a <- ggplot(df3, aes(x = est, y = Guild)) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), colour = "grey60", height = 0.175, size = 1.75) +
  geom_point(size = 4, colour = "black", fill = "grey60", shape = 21, stroke = 0.8) +
  geom_vline(xintercept = 0, linetype = 2, colour = "steelblue", size = 1) +
  scale_y_discrete(NULL) +
  labs(x = expression("Effect of bioclimatic variable on herbivory ("*z[r]*")")) +
  theme_classic(base_size = 14) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x  = element_text(size = 14, hjust = 1),
        axis.text.y  = element_text(size = 14, hjust = 1),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = "none") +
  scale_x_continuous(limits = c(-0.5, 0.5), labels = number_format(accuracy = 0.1))

p4a
ggsave("Figure 4a.pdf", p4a, height = 250, width = 125, units = "mm")

# Figure 4b
bio_vars <- paste0("Bio_", 1:19)
coef_mod2 <- lapply(bio_vars, function(bio){
  mdl  <- get(paste0("rm2_", bio), envir = .GlobalEnv)
  broom.mixed::tidy(mdl, conf.int = TRUE) %>% 
    filter(str_detect(term, "^Plant_status")) %>%
    mutate(Bio    = bio,
           Status = if_else(str_detect(term, "Native$"), "Native", "Non-native"),
           est    = estimate,
           lo     = conf.low,
           hi     = conf.high) %>%
    dplyr::select(Bio, Status, est, lo, hi)
}) %>% bind_rows()

bio_order <- rev(bio_vars)
coef_mod2 <- coef_mod2 %>%
  mutate(Bio    = factor(Bio, levels = bio_order),
         y0     = as.numeric(Bio),
         y_pos  = if_else(Status == "Native",  y0 + 0.2,
                          if_else(Status == "Non-native", y0 - 0.2, NA_real_)))

pal <- c("Native" = "#6098f9", "Non-native" = "#fb6f66")

p4b <- ggplot(coef_mod2, aes(x = est, y = y_pos, colour = Status, fill = Status)) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.175, size = 1.75) +
  geom_point(shape = 21, colour = "black", stroke = 1, size = 4) +
  geom_vline(xintercept = 0, linetype = 2, colour = "steelblue", size = 1) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values   = pal) +
  scale_y_continuous(breaks = seq_along(bio_order),
                     labels = bio_order,
                     expand = expansion(add = 0.5)) +
  labs(x = expression("Effect of bioclimatic variable on herbivory ("*z[r]*")"),
       y = NULL) +
  theme_classic(base_size = 14) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x  = element_text(size = 14, hjust = 1),
        axis.text.y  = element_text(size = 14, hjust = 1),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = "none") +
  scale_x_continuous(limits = c(-0.5, 0.5), labels = number_format(accuracy = 0.1))

p4b
ggsave("Figure 4b.pdf", p4b, height = 250, width = 125, units = "mm")
