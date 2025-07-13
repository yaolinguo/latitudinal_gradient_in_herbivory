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





### Global map ###
setwd("/Users/yaolin/Desktop/My papers/Manuscripts/2025 - Guo - BioRxiv - Meta for herbivory pattern/Submission to PNAS - 0527")
data <- read_excel("Raw_data_zr - 0705.xlsx")
data <- data[,c("Plant_status",
                "Longitude",
                "Latitude")]

data <- data %>% filter(!is.na(Longitude))
data <- data %>% filter(!is.na(Latitude))
world <- map_data("world")

p1 <- ggplot(data, aes(x = Longitude, 
                       y = Latitude, 
                       color = Plant_status, 
                       fill  = Plant_status)) +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), 
               fill = "grey", color = "grey") +
  # Native plants
  geom_point(data = data[data$Plant_status == "Native",], 
             aes(x = Longitude, y = Latitude), 
             size = 3, shape = 21, fill = "#6098f9", color = "black", 
             stroke = 0.25, alpha = 0.6) +
  # Non-native plants
  geom_point(data = data[data$Plant_status == "Non-native",], 
             aes(x = Longitude, y = Latitude), 
             size = 3, shape = 21, fill = "#fb6f66", color = "black", 
             stroke = 0.25, alpha = 0.6) + 
  xlab("Longitude (°)") + 
  ylab("Latitude (°)") +
  theme(axis.text       = element_text(size = 14),
        axis.title.x    = element_text(size = 14),
        axis.title.y    = element_text(size = 14),
        legend.title    = element_text(size = 14),
        legend.text     = element_text(size = 14),
        panel.border    = element_rect(color="black", fill=NA, size=1),
        panel.background= element_blank(),
        panel.grid.major   = element_line(colour = "grey70", linetype = "dashed", linewidth = 0.3),
        panel.grid.minor   = element_line(colour = "grey85", linetype = "dotted", linewidth = 0.2),
        legend.position = "none") +
  coord_fixed(ratio = 1.3, xlim = c(-180, 180))
p1
ggsave('./Figure 2a - 0705.pdf', p1, height = 100, width = 150, units = "mm")





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
data_zr

sample.sizes <- as.data.frame(table(data_zr$Plant_status))
sample.sizes
sample.sizes <- as.data.frame(table(data_zr$Reference))
sample.sizes

n_literatures <- length(unique(data_zr$Reference))
n_literatures
n_literatures <- length(unique(data$Reference))
n_literatures

write.csv(data_zr, file = "Raw_data_zr.csv", row.names = FALSE)

data_zr <- data_zr %>%
  mutate(zr  = as.numeric(zr),
         var = as.numeric(var))

rm1 <- rma.mv(yi     = zr,
              V      = var,
              mods   = ~ 1,
              random = ~ 1 | Plant_species + Reference,
              method = "REML",
              data   = data_zr)
rm1

rm2 <- rma.mv(yi     = zr,
              V      = var,
              mods   = ~ Plant_status - 1,
              random = ~ 1 | Plant_species + Reference,
              method = "REML",
              data   = data_zr)
rm2

sample.sizes <- as.data.frame(table(data_zr$Plant_status))
sample.sizes

coef_ci <- function(model, label){
  est  <- coef(model)
  se   <- sqrt(diag(vcov(model)))
  tibble(Type  = label,
         coef  = est,
         lower = est - 1.96 * se,
         upper = est + 1.96 * se)
}

dat <- bind_rows(coef_ci(rm1, "Total"),
                 coef_ci(rm2, c("Native", "Non-native"))) %>% 
  mutate(Type = factor(Type, levels = c("Non-native", "Native", "Total")))

data_total <- data_zr %>% 
  mutate(Plant_status = factor(Plant_status, levels = c("Non-native", "Native"))) %>% 
  bind_rows(data_zr %>% mutate(Plant_status = "Total"))

pal <- c("Total"      = "grey60",
         "Native"     = "#6098f9",
         "Non-native" = "#fb6f66")

p2 <- ggplot(dat) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, y = Type, colour = Type),
                 height = 0.075, size = 1.75) +
  geom_point(aes(x = coef, y = Type, fill = Type),
             shape = 21, colour = "black", stroke = 1, size = 4) +
  geom_vline(xintercept = 0, linetype = 2,
             colour = "steelblue", size = 1) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values   = pal) +
  scale_y_discrete(NULL) +
  labs(x = expression("Effect of latitude on herbivory ("*z[r]*")")) +
  theme_classic(base_size = 14) +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        axis.text.x = element_text(size  = 14, hjust = 1),
        axis.text.y = element_text(size  = 14, hjust = 1),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.position = "none") +
  xlab(expression("Effect of latitude on herbivory (" * z[r] * ")")) +
  scale_x_continuous(limits = c(-0.3, 0.3), labels = scales::number_format(accuracy = 0.1))
p2

ggsave("./Figure 2b - 0705.pdf", p2, height = 125, width = 100, units = "mm")
