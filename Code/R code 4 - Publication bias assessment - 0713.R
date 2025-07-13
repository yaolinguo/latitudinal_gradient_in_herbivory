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
str(data_selected$IF)
data_selected$IF <- as.numeric(data_selected$IF)

group_vars <- c("Reference", "Journal_name", "Publication_year", "IF", "Experiment_year",
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

### Test the bias ###
data <- data_zr
standard.model <- rma(zr, var, 
                      mods = ~ 1, 
                      data = data)
regtest(standard.model)
funnel(standard.model)

tf  <- trimfill(standard.model, estimator = "L0")
summary(tf)

forest.model <- rma.mv(zr, var,
                       mods = ~ 1,
                       random = list(~ 1 | Plant_species),
                       method = "REML",
                       data = data)
forest.model 

# Obtain residuals
resstandards <- rstandard.rma.mv(forest.model, type = "response")

# Obtain grand mean effect size 
grand.mean <- as.numeric(forest.model$b) 

# Create new df with residuals replacing raw
df.forest.model <- data
df.forest.model$Effect_size <- resstandards$resid + grand.mean 
df.forest.model$sei <- resstandards$se

# Funnel plot for all outcome classes
make.funnel <- function(dataset, model){
  apatheme <- theme_bw() +  
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(),
          text = element_text(family = 'Times'),
          legend.position = 'none')
  estimate <- model$b %>% as.numeric()
  SE <- model$se
  se.seq <- seq(0, max(sqrt(dataset$var)), 0.001)
  dfCI <- data.frame(ll95 = estimate - (1.96 * se.seq), 
                     ul95 = estimate + (1.96 * se.seq), 
                     ll99 = estimate - (3.29 * se.seq), 
                     ul99 = estimate + (3.29 * se.seq), 
                     se.seq = se.seq, 
                     meanll95 = estimate - (1.96 * SE), 
                     meanul95 = estimate + (1.96 * SE))
  ggplot(data, aes(x = sqrt(var), y = zr)) +
    geom_point(fill = "grey60", size = 4, shape = 21, color= "grey20", alpha = 0.75, stroke = 1) +
    xlab("Standard Error") + ylab("Effect Size (Zr)") +
    geom_line(aes(x = se.seq, y = ll95), linetype = 'dotted', data = dfCI) + # confidence lines
    geom_line(aes(x = se.seq, y = ul95), linetype = 'dotted', data = dfCI) +
    geom_line(aes(x = se.seq, y = ll99), linetype = 'dashed', data = dfCI) +
    geom_line(aes(x = se.seq, y = ul99), linetype = 'dashed', data = dfCI) +
    geom_segment(aes(x = min(se.seq), y = meanll95, xend = max(se.seq), yend = meanll95), linetype='dotdash', data=dfCI, colour = "steelblue", size =0.75) +
    geom_segment(aes(x = min(se.seq), y = meanul95, xend = max(se.seq), yend = meanul95), linetype='dotdash', data=dfCI, colour = "steelblue",size=0.75) +
    scale_x_reverse() +
    coord_flip() +
    scale_fill_brewer(palette = "Set1")+
    theme_bw() +
    theme(panel.spacing = unit(0.5, "lines"),
          panel.border= element_blank(),
          text = element_text(size = 14),
          axis.line=element_line(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.text = element_text(size = 14),
          legend.title=element_text(size = 14,),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14)) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.01))
}

funnel.plot <- make.funnel(df.forest.model, forest.model)
funnel.plot
ggsave('Figure 5a.pdf', funnel.plot, height = 100, width = 125, units = c("mm"))


# the impact by IF

fitIF <- lm(zr ~ log(IF), data = data)
summary(fitIF)

IF.plot <- ggplot(data = data, aes(x = IF, y = zr)) +
  geom_point(fill = "grey60", size = 4, shape = 21, color= "grey20", alpha = 0.75, stroke = 1) +
  geom_hline(yintercept = 0, linetype = 2, colour = "steelblue", size = 0.75) +
  scale_x_log10(limits = c(-5, 15), breaks = c(0, 1, 2, 5, 10, 15)) +
  labs(size = 'Weight (%)', y="Effect size (Zr)", x= 'Journal Impact Factor (log-scale)') + 
  guides(size = F) +
  theme_bw()+
  theme(panel.spacing = unit(0.5, "lines"),
        panel.border= element_blank(),
        text = element_text(size = 14),
        axis.line=element_line(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.text = element_text(size = 14),
        legend.title=element_text(size = 14,),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))

IF.plot
ggsave('Figure 5b.pdf', IF.plot, height = 100, width = 125, units = c("mm"))

# the impact by published year
fityear <- lm(zr ~ Publication_year, data = data)
summary(fityear)

time.plot <- data %>% 
  ggplot(aes(x = Publication_year, y= zr)) +
  geom_point(fill = "grey60", size = 4, shape = 21, color= "grey20", alpha = 0.75, stroke = 1) +
  geom_hline(yintercept = 0, linetype = 2, colour = "steelblue", size = 0.75) +
  guides(size = F) +
  labs(y="Effect size (Zr)", x= 'Year of publication') +
  theme_bw() +
  geom_smooth(method     = "lm",
              formula    = y ~ x,
              se         = TRUE,
              colour     = "grey60",
              fill       = "grey60",
              size       = 1) + 
  theme(panel.spacing      = unit(0.5, "lines"),
        panel.border       = element_blank(),
        text               = element_text(size = 14),
        axis.line=element_line(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.text        = element_text(size = 14),
        legend.title       = element_text(size = 14),
        axis.title.x       = element_text(size = 14),
        axis.title.y       = element_text(size = 14)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))
time.plot
ggsave('Figure 5c.pdf', time.plot, height = 100, width = 125,units = c("mm"))

