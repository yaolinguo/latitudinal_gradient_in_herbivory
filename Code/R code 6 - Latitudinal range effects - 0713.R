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
data_selected <- data[, selected_cols]

group_vars <- c("Reference", "Journal_name", "Publication_year", "Experiment_year",
                "Data_details", "Variable_name", "Variable_identity", "Plant_species",
                "Plant_status", "Herbivore_feeding_guild")

data_selected <- data_selected %>%
  group_by(across(all_of(group_vars))) %>%
  mutate(Latitude_extension = abs(max(Latitude,  na.rm = TRUE) - min(Latitude,  na.rm = TRUE)),
         Longitude_extension = abs(max(Longitude, na.rm = TRUE) - min(Longitude, na.rm = TRUE))) %>%
  ungroup()  

group_vars <- c("Reference", "Journal_name", "Publication_year", "Experiment_year",
                "Data_details", "Variable_name", "Variable_identity", "Plant_species",
                "Plant_status", "Herbivore_feeding_guild",
                "Latitude_extension", "Longitude_extension")
data_zr <- data_selected %>%
  filter(!is.na(Latitude) & !is.na(Herbivory)) %>%
  group_by(across(all_of(group_vars))) %>%
  summarise(cor = cor(Latitude, Herbivory, use = "complete.obs", method = "pearson"),
            n = sum(!is.na(Latitude) & !is.na(Herbivory)),
            .groups = "drop") %>%
  mutate(zr = atanh(cor),
         var = ifelse(n == 3, 1, 1 / (n - 3)))
data_zr



### Figure S2 Latitude/Longitude extension ###

# Latitude extension

data <- data_zr
data$Latitude_extension <- as.numeric(data$Latitude_extension)

data1 <- data[data$Plant_status != "Native", ]
sample.sizes.plant_identity.class <- as.data.frame(table(data1$Plant_status))
sample.sizes.plant_identity.class
fit1 <- lm(formula = zr ~ Latitude_extension - 1, data = data1)
summary(fit1)

data2 <- data[data$Plant_status != "Non-native", ]
sample.sizes.plant_identity.class <- as.data.frame(table(data2$Plant_status))
sample.sizes.plant_identity.class
fit2 <- lm(formula = zr ~ Latitude_extension - 1, data = data2)
summary(fit2)



# drawn the figure
p <- ggplot(data, aes(x = Latitude_extension, y = zr, color = Plant_status, fill = Plant_status)) +
  geom_point(data = data[data$Plant_status == "Native",], aes(x = Latitude_extension, y = zr), 
             size = 4, shape = 21, fill = "#6098f9", color = "black", stroke = 0.25, alpha = 0.7) +
  geom_point(data = data[data$Plant_status == "Non-native",], aes(x = Latitude_extension, y = zr), 
             size = 4, shape = 21, fill = "#fb6f66", color = "black", stroke = 0.25, alpha = 0.7) +  
  xlab("Latitude range (°)") +
  ylab("Herbivory effect size (Zr)") +
  theme_classic() +
  geom_smooth(data = subset(data, Plant_status == "Native"),
              method = "lm",
              formula = y ~ x-1,
              se = TRUE,
              color = "#6098f9",
              fill = "#6098f9",
              alpha = 0.2) +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        axis.text.x = element_text(size = 14, hjust = 1),
        axis.text.y = element_text(size = 14, hjust = 1),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = "none") +
  geom_hline(yintercept = 0, linetype = 2, colour = "steelblue", size = 0.75) +
  scale_x_continuous(labels = scales::number_format(accuracy = 1)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
p
ggsave('./Figure S2.pdf', p, height = 125, width = 150, units = c("mm"))