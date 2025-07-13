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
library(geodata)
library(multcompView)
library(multcomp)
library(rotl)
library(psych)
library(meta)
library(devtools)
library(processx)
library(gghalves)
library(writexl)
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
str(data_selected$IF)
data_selected$IF <- as.numeric(data_selected$IF)

group_vars <- c("Reference", "Journal_name", "Publication_year", "Experiment_year", "IF",
                "Data_details", "Variable_name", "Variable_identity", "Plant_species",
                "Plant_status", "Herbivore_feeding_guild")

data_zr <- data_selected %>%
  filter(!is.na(Latitude) & !is.na(Herbivory)) %>%
  group_by(across(all_of(group_vars))) %>%
  summarise(cor = cor(Latitude, Herbivory, use = "complete.obs", method = "pearson"),
            n = sum(!is.na(Latitude) & !is.na(Herbivory)),
            .groups = "drop") %>%
  mutate(zr = atanh(cor),
         var = ifelse(n == 3, 1, 1 / (n - 3)))

head(data_zr)
data_zr <- data_zr %>% mutate(se_zr = sqrt(var))
write.csv(data_zr, file = "data_zr.csv", row.names = FALSE)

sample.sizes.plant_identity.class <- as.data.frame(table(data_zr$Plant_status))
sample.sizes.plant_identity.class
sample.sizes.plant_identity.class <- as.data.frame(table(data_zr$Plant_status, data_zr$Herbivore_feeding_guild))
sample.sizes.plant_identity.class
sample.sizes.plant_identity.class <- as.data.frame(table(data_zr$Herbivore_feeding_guild))
sample.sizes.plant_identity.class





# REML modeling

rm3 <- rma.mv(yi     = zr,
              V      = var,
              mods   = ~ Herbivore_feeding_guild - 1,
              random = ~ 1 | Plant_species + Reference,
              method = "REML",
              data   = data_zr)
rm3

rm4 <- rma.mv(yi     = zr,
              V      = var,
              mods   = ~ paste(Herbivore_feeding_guild, Plant_status) - 1,
              random = ~ 1 | Plant_species + Reference, 
              method = "REML",
              data   = data_zr)
rm4





# Figure 3a
df3 <- broom.mixed::tidy(rm3, conf.int = TRUE) %>%
  filter(grepl("^Herbivore_feeding_guild", term)) %>%
  transmute(Guild = sub("^Herbivore_feeding_guild", "", term),
            est   = estimate,
            lo    = conf.low,
            hi    = conf.high) %>%
  mutate(Guild = factor(Guild, levels = c("Stem feeders",
                                          "Sap feeders",
                                          "Seed feeders",
                                          "Miners",
                                          "Gallers",
                                          "Grazers",
                                          "Defoliators",
                                          "All folivores"))) %>%
  arrange(Guild)

p3 <- ggplot(df3, aes(x = est, y = Guild)) +
  geom_errorbarh(aes(xmin = lo, xmax = hi),
                 colour = "grey60", height = 0.1, size = 1.75) +
  geom_point(size = 4, colour = "black", fill = "grey60", shape = 21, stroke = 0.8) +
  geom_vline(xintercept = 0, linetype = 2, colour = "steelblue", size = 1) +
  scale_y_discrete(NULL) +
  labs(x = expression("Effect of latitude on herbivory ("*z[r]*")")) +
  theme_classic(base_size = 14) +
  theme_classic() +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        axis.text.x = element_text(size  = 14, hjust = 1),
        axis.text.y = element_text(size  = 14, hjust = 1),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = "none") +
  scale_x_continuous(limits = c(-1.5, 1.5), labels = scales::number_format(accuracy = 0.1))
p3
ggsave("Figure 3a.pdf", p3, height = 150, width = 120, units = "mm")





# Figure 3b
coef4 <- broom.mixed::tidy(rm4, conf.int = TRUE) %>%
  filter(str_detect(term, "^paste\\(")) %>%
  mutate(tmp   = str_remove(term, "^paste\\([^,]+,[^\\)]+\\)"),
         Status = str_extract(tmp, "Native$|Non-native$"),
         Guild  = str_trim(str_remove(tmp, "Native$|Non-native$")),
         est    = estimate,
         lo     = conf.low,
         hi     = conf.high) %>%
  dplyr::select(Guild, Status, est, lo, hi)
guild_order <- c("All folivores",
                 "Defoliators",
                 "Grazers",
                 "Gallers",
                 "Miners",
                 "Seed feeders",
                 "Sap feeders",
                 "Stem feeders")
coef4 <- coef4 %>%
  mutate(Guild  = factor(Guild,  levels = guild_order),
         Status = factor(Status, levels = c("Native","Non-native")),
         y0    = as.numeric(Guild),
         y_pos = if_else(Status == "Native",  y0 - 0.2, 
                         if_else(Status == "Non-native", y0 + 0.2, NA_real_)))

pal <- c("Native"     = "#6098f9",
         "Non-native" = "#fb6f66")

p4 <- ggplot(coef4, aes(x=est, y=y_pos, colour=Status, fill=Status)) +
  geom_errorbarh(aes(xmin=lo, xmax=hi), height = 0.1, size = 1.75) +
  geom_point(shape = 21, colour = "black", stroke = 1, size = 4) +
  geom_vline(xintercept=0, linetype="dashed", colour="steelblue") +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values   = pal) +
  scale_y_continuous(breaks = seq_along(guild_order),
                     labels = guild_order,
                     expand = expansion(add = 0.5),
                     trans = "reverse") +
  labs(x = expression("Effect of latitude on herbivory ("*z[r]*")"), y = NULL) +
  theme_classic(base_size = 14) +
  theme_classic() +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        axis.text.x = element_text(size  = 14, hjust = 1),
        axis.text.y = element_text(size  = 14, hjust = 1),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = "none") +
  scale_x_continuous(limits = c(-1.5, 1.5), labels = scales::number_format(accuracy = 0.1))
p4

ggsave("Figure 4b.pdf", p4, height = 150, width = 120, units = "mm")

