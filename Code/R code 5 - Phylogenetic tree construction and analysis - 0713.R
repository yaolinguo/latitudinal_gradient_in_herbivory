#########################################################################################################
### Guo et al., 2025. Herbivory increases towards lower latitudes in native but not introduced plants ###
#########################################################################################################

library(taxize)
library(ggtree)
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
library(V.PhyloMaker)

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
  group_by(across(all_of(group_vars))) %>%
  summarise(cor = cor(Latitude, Herbivory, use = "complete.obs", method = "pearson"),
            n = sum(!is.na(Latitude) & !is.na(Herbivory)),
            .groups = "drop") %>%
  mutate(zr = atanh(cor),
         var = ifelse(n == 3, 1, 1 / (n - 3)))





### Tree and pholy relatinship ###
data <- data_zr
number_of_species <- data %>% distinct(Plant_species) %>% nrow()
number_of_species

data <- data[,c("Reference",
                "Plant_status",
                "Plant_species",
                "zr")]
table(data$Plant_status)

sample.sizes.plant_identity.class <- as.data.frame(table(data$Plant_species))
sample.sizes.plant_identity.class

data %>% distinct(Plant_status, Plant_species) %>%
  group_by(Plant_status) %>%
  summarise(num_species = n_distinct(Plant_species))

overlaps <- data %>%
  distinct(Plant_species, Plant_status) %>%
  group_by(Plant_species) %>%
  summarise(num_identity = n_distinct(Plant_status),
            identities   = paste(unique(Plant_status), collapse = ", ")) %>%
  filter(num_identity > 1)
overlaps

data <- data %>%  mutate(Plant_species = str_trim(Plant_species))
synonyms <- c("Betula pubscens"           = "Betula pubescens",
              "Xylocarpus australiasicus" = "Xylocarpus australasicus",
              "Fallopia japonica"         = "Reynoutria japonica",
              "Padus avium"               = "Prunus padus",
              "Sapium sebiferum"          = "Triadica sebifera")
data <- data %>% mutate(Plant_species = ifelse(Plant_species %in% names(synonyms),
                                               synonyms[Plant_species],
                                               Plant_species))
unique(data$Plant_species)
length(unique(data$Plant_species))
length(unique(data$Reference))

df <- data %>%  distinct(Plant_species) %>% arrange(Plant_species) %>% rename(Species = Plant_species)
print(df, n = 140)

Family <- c(
  "Fabaceae",        # 1  Acacia falcata
  "Sapindaceae",     # 2  Acer platanoides
  "Sapindaceae",     # 3  Acer rubrum
  "Sapindaceae",     # 4  Acer saccharum
  "Asteraceae",      # 5  Achillea millefolium
  "Primulaceae",     # 6  Aegiceras corniculatum
  "Betulaceae",      # 7  Alnus incana
  "Betulaceae",      # 8  Alnus viridis
  "Amaranthaceae",   # 9  Alternanthera philoxeroides
  "Amaranthaceae",   # 10 Alternanthera sessilis
  "Asteraceae",      # 11 Ambrosia artemisiifolia
  "Asteraceae",      # 12 Anaphalis margaritacea
  "Poaceae",         # 13 Andropogon gerardii
  "Myrtaceae",       # 14 Angophora floribunda
  "Myrtaceae",       # 15 Angophora hispida
  "Brassicaceae",    # 16 Arabidopsis lyrata
  "Asteraceae",      # 17 Arctium minus
  "Primulaceae",     # 18 Ardisia elliptica
  "Primulaceae",     # 19 Ardisia escallonioides
  "Acanthaceae",     # 20 Avicennia marina
  "Berberidaceae",   # 21 Berberis thunbergii
  "Betulaceae",      # 22 Betula glandulosa
  "Betulaceae",      # 23 Betula nana
  "Betulaceae",      # 24 Betula pendula
  "Betulaceae",      # 25 Betula platyphylla
  "Betulaceae",      # 26 Betula pubescens
  "Poaceae",         # 27 Bouteloua dactyloides
  "Poaceae",         # 28 Bouteloua eriopoda
  "Poaceae",         # 29 Bouteloua gracilis
  "Rhizophoraceae",  # 30 Bruguiera gymnorrhiza
  "Theaceae",        # 31 Camellia japonica
  "Campanulaceae",   # 32 Campanula americana
  "Betulaceae",      # 33 Carpinus betulus
  "Celastraceae",    # 34 Celastrus orbiculatus
  "Asteraceae",      # 35 Cirsium arvense
  "Betulaceae",      # 36 Corylus avellana
  "Rosaceae",        # 37 Cotoneaster melanocarpus
  "Thymelaeaceae",   # 38 Daphne laureola
  "Rosaceae",        # 39 Dasiphora fruticosa
  "Fabaceae",        # 40 Daviesia corymbosa
  "Atherospermataceae", # 41 Doryphora sassafras
  "Elaeagnaceae",    # 42 Elaeagnus umbellata
  "Myrtaceae",       # 43 Eucalyptus camaldulensis
  "Myrtaceae",       # 44 Eucalyptus globulus
  "Myrtaceae",       # 45 Eucalyptus melliodora
  "Myrtaceae",       # 46 Eucalyptus populnea
  "Asteraceae",      # 47 Eupatorium maculatum
  "Euphorbiaceae",   # 48 Excoecaria agallocha
  "Fagaceae",        # 49 Fagus crenata
  "Fagaceae",        # 50 Fagus grandifolia
  "Fagaceae",        # 51 Fagus sylvatica
  "Malvaceae",       # 52 Gossypium hirsutum
  "Ranunculaceae",   # 53 Helleborus foetidus
  "Aquifoliaceae",   # 54 Ilex aquifolium
  "Cupressaceae",    # 55 Juniperus communis
  "Brassicaceae",    # 56 Lepidium draba
  "Asteraceae",      # 57 Leucanthemum vulgare
  "Altingiaceae",    # 58 Liquidambar styraciflua
  "Lythraceae",      # 59 Lythrum salicaria
  "Asteraceae",      # 60 Mikania micrantha
  "Nothofagaceae",   # 61 Nothofagus pumilio
  "Onagraceae",      # 62 Oenothera biennis
  "Tetrameristaceae",# 63 Pelliciera rhizophorae
  "Poaceae",         # 64 Phragmites australis
  "Phytolaccaceae",  # 65 Phytolacca americana
  "Phytolaccaceae",  # 66 Phytolacca rivinoides
  "Pinaceae",        # 67 Picea abies
  "Pinaceae",        # 68 Pinus sylvestris
  "Piperaceae",      # 69 Piper aduncum
  "Piperaceae",      # 70 Piper aequale
  "Polygonaceae",    # 71 Polygonum macrophyllum
  "Salicaceae",      # 72 Populus tremula
  "Proteaceae",      # 73 Protea laurifolia
  "Proteaceae",      # 74 Protea magnifica
  "Rosaceae",        # 75 Prunus padus
  "Rosaceae",        # 76 Prunus serotina
  "Fagaceae",        # 77 Quercus acutissima
  "Fagaceae",        # 78 Quercus alba
  "Fagaceae",        # 79 Quercus aliena
  "Fagaceae",        # 80 Quercus bicolor
  "Fagaceae",        # 81 Quercus castaneifolia
  "Fagaceae",        # 82 Quercus cerris
  "Fagaceae",        # 83 Quercus dentata
  "Fagaceae",        # 84 Quercus gambelii
  "Fagaceae",        # 85 Quercus garryana
  "Fagaceae",        # 86 Quercus georgiana
  "Fagaceae",        # 87 Quercus glandulifera
  "Fagaceae",        # 88 Quercus grisea
  "Fagaceae",        # 89 Quercus ilicifolia
  "Fagaceae",        # 90 Quercus imbricaria
  "Fagaceae",        # 91 Quercus lyrata
  "Fagaceae",        # 92 Quercus michauxii
  "Fagaceae",        # 93 Quercus mongolica
  "Fagaceae",        # 94 Quercus muehlenbergii
  "Fagaceae",        # 95 Quercus phellos
  "Fagaceae",        # 96 Quercus prinoides
  "Fagaceae",        # 97 Quercus pubescens
  "Fagaceae",        # 98 Quercus robur
  "Fagaceae",        # 99 Quercus rubra
  "Fagaceae",        # 100 Quercus rugosa
  "Fagaceae",        # 101 Quercus serrata
  "Fagaceae",        # 102 Quercus shumardii
  "Fagaceae",        # 103 Quercus suber
  "Fagaceae",        # 104 Quercus texana
  "Fagaceae",        # 105 Quercus turbinella
  "Fagaceae",        # 106 Quercus variabilis
  "Fagaceae",        # 107 Quercus velutina
  "Fagaceae",        # 108 Quercus virginiana
  "Rhizophoraceae",  # 109 Rhizophora mangle
  "Rhizophoraceae",  # 110 Rhizophora stylosa
  "Asteraceae",      # 111 Rudbeckia hirta
  "Acanthaceae",     # 112 Ruellia nudiflora
  "Alismataceae",    # 113 Sagittaria latifolia
  "Salicaceae",      # 114 Salix caprea
  "Salicaceae",      # 115 Salix dasyclados
  "Salicaceae",      # 116 Salix lanata
  "Salicaceae",      # 117 Salix myrtilloides
  "Salicaceae",      # 118 Salix phylicifolia
  "Salicaceae",      # 119 Salix polaris
  "Asteraceae",      # 120 Saussurea pulchra
  "Poaceae",         # 121 Schizachyrium scoparium
  "Asteraceae",      # 122 Senecio madagascariensis
  "Asteraceae",      # 123 Senecio pinnatifolius
  "Solanaceae",      # 124 Solanum carolinense
  "Asteraceae",      # 125 Solidago altissima
  "Asteraceae",      # 126 Solidago canadensis
  "Asteraceae",      # 127 Sonchus arvensis
  "Rosaceae",        # 128 Sorbus aucuparia
  "Poaceae",         # 129 Spartina alterniflora
  "Asteraceae",      # 130 Taraxacum mongolicum
  "Asteraceae",      # 131 Taraxacum officinale
  "Malvaceae",       # 132 Tilia cordata
  "Euphorbiaceae",   # 133 Triadica sebifera
  "Ulmaceae",        # 134 Ulmus glabra
  "Ericaceae"        # 135 Vaccinium myrtillus
)

df <- df %>% mutate(Family = Family)
print(df, n = 140)

data <- data %>%
  group_by(Plant_species) %>%
  mutate(Plant_status = case_when(n_distinct(Plant_status) == 2 ~ "Overlap",
                                  Plant_status == "Native"     ~ "Native",
                                  Plant_status == "Non-native" ~ "Non-native")) %>%
  ungroup()

data <- data %>%
  group_by(Plant_species) %>%
  mutate(zr = mean(zr, na.rm = TRUE)) %>%
  slice(1) %>%
  ungroup()
table(data$Plant_status)

data <- data %>% left_join(df, by = c("Plant_species" = "Species"))

data_tree <- data[, c("zr", "Plant_species", "Family", "Plant_status")] %>%
  mutate(Plant_species = str_replace_all(Plant_species, "﻿", "")) %>%
  separate(Plant_species, into = c("genus", "species"), sep = " ", remove = FALSE) %>%
  filter(species != "spp.")

data_tree <- data_tree %>% mutate(species = if_else(Plant_species == "Quercus alba",
                                                    "Quercus alba",
                                                    species))
data_tree <- data_tree %>% mutate(species = if_else(Plant_species == "Campanula americana",
                                                    "Campanula americana",
                                                    species))
data_tree <- data_tree %>% mutate(species = if_else(Plant_species == "Camellia japonica",
                                                    "Camellia japonica",
                                                    species))
data_tree <- data_tree %>% mutate(species = if_else(Plant_species == "Arabidopsis lyrata",
                                                    "Arabidopsis lyrata",
                                                    species))
data_tree <- data_tree %>% mutate(species = if_else(Plant_species == "Betula pubescens",
                                                    "Betula pubescens",
                                                    species))
  
species.list <- data_tree %>% dplyr::select(species, genus, Family)
phylo_tree_result <- phylo.maker(sp.list = species.list)
phylo_tree <- phylo_tree_result$scenario.3
tree_data <- phylo_tree
tree_data$tip.label <- gsub("_", " ", tree_data$tip.label)
data_tree <- data_tree %>% rename(label = species)
data_tree <- data_tree[ , !duplicated(names(data_tree))]
unique(tree_data$tip.label)
unique(data_tree$label)

p <- ggtree(tree_data, layout = "circular") + geom_tree()

p$data$label <- tolower(p$data$label)
data_tree$label <- tolower(data_tree$label)
p$data <- p$data %>% left_join(data_tree, by = "label")
p$data <- p$data %>% mutate(species_short = word(label, -1)) 
p$data$species_short <- str_to_title(p$data$species_short)
filtered_data <- p$data %>% filter(species_short %in% tree_data$tip.label)

p <- p %<+% filtered_data
p1 <- p +
  geom_tippoint(aes(color = Plant_status, size = zr)) +
  geom_tiplab(aes(label = Plant_species, color = Plant_status),
              size      = 3,
              align     = TRUE,
              offset    = 25,
              linetype  = "dotted") +
  theme(plot.margin = unit(c(3, 3, 3, 3), "cm"),
        text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = c(1.9, 0.2),
        legend.justification = c(1.9, 0.2)) +
  scale_color_manual(values = c("Native" = "#6098f9", "Non-native" = "#fb6f66")) +
  scale_size_continuous(range = c(-0.3, 2))
p1

ggsave('./Figure S1.pdf', p1, height = 200, width = 500, units = c("mm"))





# test for overall effect sizes
matched_list <- phylo_tree_result$species.list
matched_list <- matched_list %>%
  mutate(species = if_else(species == "Quercus alba", "alba", species)) %>%
  mutate(species = if_else(species == "Betula pubescens", "pubescens", species)) %>%
  mutate(species = if_else(species == "Camellia japonica", "japonica", species)) %>%
  mutate(species = if_else(species == "Campanula americana", "americana", species)) %>%
  mutate(species = if_else(species == "Arabidopsis lyrata", "lyrata", species))
genus_species_names <- paste0(matched_list$genus, "_", matched_list$species)
phylo_tree$tip.label <- genus_species_names
data$Plant_species <- gsub(" ", "_", data$Plant_species)
common <- intersect(phylo_tree$tip.label, data$Plant_species)
cat("Matched species count:", length(common), "\n")
print(common)

if (is.null(phylo_tree$edge.length)) {
  cat("Tree has no branch lengths. Setting them all to 1.\n")
  phylo_tree$edge.length <- rep(1, nrow(phylo_tree$edge))
}

effect_size_species <- data %>%
  dplyr::group_by(Plant_species) %>%
  dplyr::summarize(zr = mean(zr, na.rm = TRUE))

effect_size_species <- effect_size_species[!is.na(effect_size_species$zr), ]
rownames(effect_size_species) <- effect_size_species$Plant_species
effect_size_species$Plant_species <- gsub(" ", "_", effect_size_species$Plant_species)
rownames(effect_size_species)     <- gsub(" ", "_", rownames(effect_size_species))

cat("Tree tip labels example:\n")
print(head(phylo_tree$tip.label, 1000))
cat("Data rownames example:\n")
print(head(rownames(effect_size_species), 10))
common_species <- intersect(phylo_tree$tip.label, rownames(effect_size_species))
cat("Number of common species between tree and data:", length(common_species), "\n")
effect_size_species <- effect_size_species[common_species, , drop = FALSE]
effect_size_species <- effect_size_species[match(phylo_tree$tip.label, rownames(effect_size_species)), ]
effect_size_species <- effect_size_species[!is.na(rownames(effect_size_species)), ]
names(effect_size_species$zr) <- rownames(effect_size_species)

lambda_result_overall <- phylosig(tree   = phylo_tree,
                              x          = effect_size_species$zr,
                              method     = "lambda",
                              test       = TRUE)
print(lambda_result_overall)





# For native species
setwd("/Users/yaolin/Desktop/My papers/Manuscripts/Guo et al., 2024 - BioRxiv - Meta/New version - 20250121")
data <- read_excel("RawDataNewVersion_field.xlsx")
write.csv(data, file = "data.csv", row.names = FALSE)
data <- read.csv("data.csv")

data<- subset(data, Plant_identity == "Native")

number_of_species <- data %>%
  distinct(Plant_species) %>%
  nrow()
number_of_species

data <- data %>% filter(!is.na(Plant_identity))
data <- data %>% filter(!is.na(Plant_species))
data <- data %>% filter(!is.na(Variable_identity))
data <- data %>% filter(!is.na(Effect_size))
data <- data[,c("Reference",
                "Plant_identity",
                "Plant_species",
                "Effect_size")]
table(data$Plant_identity)

sample.sizes.plant_identity.class <- as.data.frame(table(data$Plant_species))
sample.sizes.plant_identity.class

data %>%
  distinct(Plant_identity, Plant_species) %>%
  group_by(Plant_identity) %>%
  summarise(num_species = n_distinct(Plant_species))
table(data$Plant_identity)

overlaps <- data %>%
  distinct(Plant_species, Plant_identity) %>%
  group_by(Plant_species) %>%
  summarise(
    num_identity = n_distinct(Plant_identity),
    identities   = paste(unique(Plant_identity), collapse = ", ")
  ) %>%
  filter(num_identity > 1)
overlaps

data <- data %>%  mutate(Plant_species = str_trim(Plant_species))

synonyms <- c(
  "Betula pubscens"           = "Betula pubescens",
  "Xylocarpus australiasicus" = "Xylocarpus australasicus",
  "Fallopia japonica" = "Reynoutria japonica",
  "Padus avium" = "Prunus padus",
  "Sapium sebiferum"          = "Triadica sebifera")
data <- data %>%
  mutate(Plant_species = ifelse(
    Plant_species %in% names(synonyms),
    synonyms[Plant_species],
    Plant_species))

unique(data$Plant_species)
length(unique(data$Plant_species))

Species <- c(
  "Acacia falcata",
  "Acacia obtusata",
  "Acanthus ilicifolius",
  "Acer platanoides",
  "Acer rubrum",
  "Acer saccharum",
  "Achillea millefolium",
  "Acrostichum speciosum",
  "Aegiceras corniculatum",
  "Aextoxicon punctatum",
  "Alnus incana",
  "Alnus viridis",
  "Alternanthera philoxeroides",
  "Alternanthera sessilis",
  "Ambrosia artemisiifolia",
  "Anaphalis margaritacea",
  "Andropogon gerardii",
  "Angophora floribunda",
  "Angophora hispida",
  "Arabidopsis lyrata",
  "Arctium minus",
  "Ardisia elliptica",
  "Ardisia escallonioides",
  "Avena sativa",
  "Avicennia germinans",
  "Avicennia marina",
  "Berberis thunbergii",
  "Betonica officinalis",
  "Betula glandulosa",
  "Betula nana",
  "Betula pendula",
  "Betula platyphylla",
  "Betula pubescens",
  "Bouteloua dactyloides",
  "Bouteloua eriopoda",
  "Bouteloua gracilis",
  "Bruguiera gymnorrhiza",
  "Camellia japonica",
  "Campanula americana",
  "Carpinus betulus",
  "Celastrus orbiculatus",
  "Chamaecrista fasciculata",
  "Cirsium arvense",
  "Coffea arabica",
  "Cordia alliodora",
  "Cornus mas",
  "Corylus avellana",
  "Cotoneaster melanocarpus",
  "Daphne laureola",
  "Dasiphora fruticosa",
  "Daviesia corymbosa",
  "Doryphora sassafras",
  "Elaeagnus umbellata",
  "Eucalyptus camaldulensis",
  "Eucalyptus globulus",
  "Eucalyptus melliodora",
  "Eucalyptus populnea",
  "Eupatorium maculatum",
  "Excoecaria agallocha",
  "Fagus crenata",
  "Fagus grandifolia",
  "Fagus sylvatica",
  "Gossypium hirsutum",
  "Gossypium thurberi",
  "Helleborus foetidus",
  "Heritiera littoralis",
  "Ilex aquifolium",
  "Juniperus communis",
  "Kandelia obovata",
  "Lepidium draba",
  "Leucanthemum vulgare",
  "Liquidambar styraciflua",
  "Lumnitzera racemosa",
  "Lythrum salicaria",
  "Mikania micrantha",
  "Nothofagus pumilio",
  "Oenothera biennis",
  "Pelliciera rhizophorae",
  "Phragmites australis",
  "Phytolacca americana",
  "Phytolacca rivinoides",
  "Picea abies",
  "Pinus sylvestris",
  "Piper aduncum",
  "Piper aequale",
  "Polygonum macrophyllum",
  "Populus tremula",
  "Protea laurifolia",
  "Protea magnifica",
  "Prunus padus",
  "Quercus acutissima",
  "Quercus alba",
  "Quercus aliena",
  "Quercus bicolor",
  "Quercus castaneifolia",
  "Quercus cerris",
  "Quercus dentata",
  "Quercus gambelii",
  "Quercus garryana",
  "Quercus georgiana",
  "Quercus glandulifera",
  "Quercus grisea",
  "Quercus ilicifolia",
  "Quercus imbricaria",
  "Quercus lyrata",
  "Quercus michauxii",
  "Quercus mongolica",
  "Quercus muehlenbergii",
  "Quercus phellos",
  "Quercus prinoides",
  "Quercus pubescens",
  "Quercus robur",
  "Quercus rubra",
  "Quercus rugosa",
  "Quercus serrata",
  "Quercus shumardii",
  "Quercus suber",
  "Quercus texana",
  "Quercus turbinella",
  "Quercus variabilis",
  "Quercus velutina",
  "Quercus virginiana",
  "Reynoutria japonica",
  "Rhizophora apiculata",
  "Rhizophora mangle",
  "Rhizophora stylosa",
  "Rudbeckia hirta",
  "Ruellia nudiflora",
  "Sagittaria latifolia",
  "Salix caprea",
  "Salix dasyclados",
  "Salix myrsinifolia",
  "Salix phylicifolia",
  "Salix polaris",
  "Saussurea pulchra",
  "Schizachyrium scoparium",
  "Senecio madagascariensis",
  "Senecio pinnatifolius",
  "Solanum carolinense",
  "Solanum dulcamara",
  "Solanum lycocarpum",
  "Solidago altissima",
  "Solidago canadensis",
  "Sonchus arvensis",
  "Sonneratia alba",
  "Sonneratia apetala",
  "Sorbus aucuparia",
  "Taraxacum mongolicum",
  "Taraxacum officinale",
  "Taxus baccata",
  "Tilia cordata",
  "Triadica sebifera",
  "Turnera ulmifolia",
  "Ulmus glabra",
  "Viburnum dilatatum",
  "Vincetoxicum nigrum",
  "Xylocarpus australasicus",
  "Xylocarpus granatum"
)

Family <- c(
  "Fabaceae",           # 1. Acacia falcata
  "Fabaceae",           # 2. Acacia obtusata
  "Acanthaceae",        # 3. Acanthus ilicifolius
  "Sapindaceae",        # 4. Acer platanoides
  "Sapindaceae",        # 5. Acer rubrum
  "Sapindaceae",        # 6. Acer saccharum
  "Asteraceae",         # 7. Achillea millefolium
  "Pteridaceae",        # 8. Acrostichum speciosum
  "Primulaceae",        # 9. Aegiceras corniculatum
  "Aextoxicaceae",      # 10. Aextoxicon punctatum
  "Betulaceae",         # 11. Alnus incana
  "Betulaceae",         # 12. Alnus viridis
  "Amaranthaceae",      # 13. Alternanthera philoxeroides
  "Amaranthaceae",      # 14. Alternanthera sessilis
  "Asteraceae",         # 15. Ambrosia artemisiifolia
  "Asteraceae",         # 16. Anaphalis margaritacea
  "Poaceae",            # 17. Andropogon gerardii
  "Myrtaceae",          # 18. Angophora floribunda
  "Myrtaceae",          # 19. Angophora hispida
  "Brassicaceae",       # 20. Arabidopsis lyrata
  "Asteraceae",         # 21. Arctium minus
  "Primulaceae",        # 22. Ardisia elliptica
  "Primulaceae",        # 23. Ardisia escallonioides
  "Poaceae",            # 24. Avena sativa
  "Acanthaceae",        # 25. Avicennia germinans
  "Acanthaceae",        # 26. Avicennia marina
  "Berberidaceae",      # 27. Berberis thunbergii
  "Lamiaceae",          # 28. Betonica officinalis
  "Betulaceae",         # 29. Betula glandulosa
  "Betulaceae",         # 30. Betula nana
  "Betulaceae",         # 31. Betula pendula
  "Betulaceae",         # 32. Betula platyphylla
  "Betulaceae",         # 33. Betula pubescens
  "Poaceae",            # 34. Bouteloua dactyloides
  "Poaceae",            # 35. Bouteloua eriopoda
  "Poaceae",            # 36. Bouteloua gracilis
  "Rhizophoraceae",     # 37. Bruguiera gymnorrhiza
  "Theaceae",           # 38. Camellia japonica
  "Campanulaceae",      # 39. Campanula americana
  "Betulaceae",         # 40. Carpinus betulus
  "Celastraceae",       # 41. Celastrus orbiculatus
  "Fabaceae",           # 42. Chamaecrista fasciculata
  "Asteraceae",         # 43. Cirsium arvense
  "Rubiaceae",          # 44. Coffea arabica
  "Boraginaceae",       # 45. Cordia alliodora
  "Cornaceae",          # 46. Cornus mas
  "Betulaceae",         # 47. Corylus avellana
  "Rosaceae",           # 48. Cotoneaster melanocarpus
  "Thymelaeaceae",      # 49. Daphne laureola
  "Rosaceae",           # 50. Dasiphora fruticosa
  "Fabaceae",           # 51. Daviesia corymbosa
  "Atherospermataceae", # 52. Doryphora sassafras
  "Elaeagnaceae",       # 53. Elaeagnus umbellata
  "Myrtaceae",          # 54. Eucalyptus camaldulensis
  "Myrtaceae",          # 55. Eucalyptus globulus
  "Myrtaceae",          # 56. Eucalyptus melliodora
  "Myrtaceae",          # 57. Eucalyptus populnea
  "Asteraceae",         # 58. Eupatorium maculatum
  "Euphorbiaceae",      # 59. Excoecaria agallocha
  "Fagaceae",           # 60. Fagus crenata
  "Fagaceae",           # 61. Fagus grandifolia
  "Fagaceae",           # 62. Fagus sylvatica
  "Malvaceae",          # 63. Gossypium hirsutum
  "Malvaceae",          # 64. Gossypium thurberi
  "Ranunculaceae",      # 65. Helleborus foetidus
  "Malvaceae",          # 66. Heritiera littoralis
  "Aquifoliaceae",      # 67. Ilex aquifolium
  "Cupressaceae",       # 68. Juniperus communis
  "Rhizophoraceae",     # 69. Kandelia obovata
  "Brassicaceae",       # 70. Lepidium draba
  "Asteraceae",         # 71. Leucanthemum vulgare
  "Altingiaceae",       # 72. Liquidambar styraciflua
  "Combretaceae",       # 73. Lumnitzera racemosa
  "Lythraceae",         # 74. Lythrum salicaria
  "Asteraceae",         # 75. Mikania micrantha
  "Nothofagaceae",      # 76. Nothofagus pumilio
  "Onagraceae",         # 77. Oenothera biennis
  "Tetrameristaceae",   # 78. Pelliciera rhizophorae
  "Poaceae",            # 79. Phragmites australis
  "Phytolaccaceae",     # 80. Phytolacca americana
  "Phytolaccaceae",     # 81. Phytolacca rivinoides
  "Pinaceae",           # 82. Picea abies
  "Pinaceae",           # 83. Pinus sylvestris
  "Piperaceae",         # 84. Piper aduncum
  "Piperaceae",         # 85. Piper aequale
  "Polygonaceae",       # 86. Polygonum macrophyllum
  "Salicaceae",         # 87. Populus tremula
  "Proteaceae",         # 88. Protea laurifolia
  "Proteaceae",         # 89. Protea magnifica
  "Rosaceae",           # 90. Prunus padus
  "Fagaceae",           # 91. Quercus acutissima
  "Fagaceae",           # 92. Quercus alba
  "Fagaceae",           # 93. Quercus aliena
  "Fagaceae",           # 94. Quercus bicolor
  "Fagaceae",           # 95. Quercus castaneifolia
  "Fagaceae",           # 96. Quercus cerris
  "Fagaceae",           # 97. Quercus dentata
  "Fagaceae",           # 98. Quercus gambelii
  "Fagaceae",           # 99. Quercus garryana
  "Fagaceae",           # 100. Quercus georgiana
  "Fagaceae",           # 101. Quercus glandulifera
  "Fagaceae",           # 102. Quercus grisea
  "Fagaceae",           # 103. Quercus ilicifolia
  "Fagaceae",           # 104. Quercus imbricaria
  "Fagaceae",           # 105. Quercus lyrata
  "Fagaceae",           # 106. Quercus michauxii
  "Fagaceae",           # 107. Quercus mongolica
  "Fagaceae",           # 108. Quercus muehlenbergii
  "Fagaceae",           # 109. Quercus phellos
  "Fagaceae",           # 110. Quercus prinoides
  "Fagaceae",           # 111. Quercus pubescens
  "Fagaceae",           # 112. Quercus robur
  "Fagaceae",           # 113. Quercus rubra
  "Fagaceae",           # 114. Quercus rugosa
  "Fagaceae",           # 115. Quercus serrata
  "Fagaceae",           # 116. Quercus shumardii
  "Fagaceae",           # 117. Quercus suber
  "Fagaceae",           # 118. Quercus texana
  "Fagaceae",           # 119. Quercus turbinella
  "Fagaceae",           # 120. Quercus variabilis
  "Fagaceae",           # 121. Quercus velutina
  "Fagaceae",           # 122. Quercus virginiana
  "Polygonaceae",       # 123. Reynoutria japonica
  "Rhizophoraceae",     # 124. Rhizophora apiculata
  "Rhizophoraceae",     # 125. Rhizophora mangle
  "Rhizophoraceae",     # 126. Rhizophora stylosa
  "Asteraceae",         # 127. Rudbeckia hirta
  "Acanthaceae",        # 128. Ruellia nudiflora
  "Alismataceae",       # 129. Sagittaria latifolia
  "Salicaceae",         # 130. Salix caprea
  "Salicaceae",         # 131. Salix dasyclados
  "Salicaceae",         # 132. Salix myrsinifolia
  "Salicaceae",         # 133. Salix phylicifolia
  "Salicaceae",         # 134. Salix polaris
  "Asteraceae",         # 135. Saussurea pulchra
  "Poaceae",            # 136. Schizachyrium scoparium
  "Asteraceae",         # 137. Senecio madagascariensis
  "Asteraceae",         # 138. Senecio pinnatifolius
  "Solanaceae",         # 139. Solanum carolinense
  "Solanaceae",         # 140. Solanum dulcamara
  "Solanaceae",         # 141. Solanum lycocarpum
  "Asteraceae",         # 142. Solidago altissima
  "Asteraceae",         # 143. Solidago canadensis
  "Asteraceae",         # 144. Sonchus arvensis
  "Lythraceae",         # 145. Sonneratia alba
  "Lythraceae",         # 146. Sonneratia apetala
  "Rosaceae",           # 147. Sorbus aucuparia
  "Asteraceae",         # 148. Taraxacum mongolicum
  "Asteraceae",         # 149. Taraxacum officinale
  "Taxaceae",           # 150. Taxus baccata
  "Malvaceae",          # 151. Tilia cordata
  "Euphorbiaceae",      # 152. Triadica sebifera
  "Passifloraceae",     # 153. Turnera ulmifolia
  "Ulmaceae",           # 154. Ulmus glabra
  "Adoxaceae",          # 155. Viburnum dilatatum
  "Apocynaceae",        # 156. Vincetoxicum nigrum
  "Meliaceae",          # 157. Xylocarpus australasicus
  "Meliaceae"           # 158. Xylocarpus granatum
)

df <- data.frame(Species, Family)
dim(df)

data <- data %>%
  left_join(
    df,
    by = c("Plant_species" = "Species")
  )

data <- data %>%
  group_by(Plant_species) %>%
  mutate(zr = mean(zr, na.rm = TRUE)) %>%
  slice(1) %>%
  ungroup()
table(data$Plant_status)

data_tree <- data[, c("zr", "Plant_species", "Family", "Plant_status")] %>%
  mutate(Plant_species = str_replace_all(Plant_species, "﻿", "")) %>%
  separate(Plant_species, into = c("genus", "species"), sep = " ", remove = FALSE) %>%
  filter(species != "spp.")
data_tree <- data_tree %>%
  mutate(
    species = if_else(
      Plant_species == "Quercus alba",
      "Quercus alba",
      species
    )
  )
data_tree <- data_tree %>%
  mutate(
    species = if_else(
      Plant_species == "Campanula americana",
      "Campanula americana",
      species
    )
  )
data_tree <- data_tree %>%
  mutate(
    species = if_else(
      Plant_species == "Camellia japonica",
      "Camellia japonica",
      species
    )
  )
data_tree <- data_tree %>%
  mutate(
    species = if_else(
      Plant_species == "Arabidopsis lyrata",
      "Arabidopsis lyrata",
      species
    )
  )
data_tree <- data_tree %>%
  mutate(
    species = if_else(
      Plant_species == "Betula pubescens",
      "Betula pubescens",
      species
    )
  )

species.list <- data_tree %>% dplyr::select(species, genus, Family)
phylo_tree_result <- phylo.maker(sp.list = species.list)
phylo_tree <- phylo_tree_result$scenario.3
tree_data <- phylo_tree
tree_data$tip.label <- gsub("_", " ", tree_data$tip.label)
data_tree <- data_tree %>% rename(label = species)
data_tree <- data_tree[ , !duplicated(names(data_tree))]
unique(tree_data$tip.label)
unique(data_tree$label)

p <- ggtree(tree_data, layout = "circular") + geom_tree()
p$data$label <- tolower(p$data$label)
data_tree$label <- tolower(data_tree$label)
p$data <- p$data %>% left_join(data_tree, by = "label")
p$data <- p$data %>% mutate(species_short = word(label, -1)) 
p$data$species_short <- str_to_title(p$data$species_short)
filtered_data <- p$data %>% filter(species_short %in% tree_data$tip.label)
p <- p %<+% filtered_data
p





# Phylogenetic distance analysis

synonyms <- c("Betula pubscens"           = "Betula pubescens",
              "Xylocarpus australiasicus" = "Xylocarpus australasicus",
              "Fallopia japonica"         = "Reynoutria japonica",
              "Padus avium"               = "Prunus padus",
              "Sapium sebiferum"          = "Triadica sebifera")

normalize_name <- function(x) gsub(" ", "_", trimws(tolower(x)))

calc_lambda <- function(dat_sub, tree_raw) {
  
  tree <- tree_raw
  tree$tip.label        <- normalize_name(tree$tip.label)
  dat_sub$Plant_species <- normalize_name(dat_sub$Plant_species)
  
  common_sp <- intersect(tree$tip.label, dat_sub$Plant_species)
  cat(" → 物种交集数:", length(common_sp), "\n")
  stopifnot(length(common_sp) >= 3)                    # 少于 3 直接报错
  tree <- drop.tip(tree, setdiff(tree$tip.label, common_sp))
  
  zr_vec <- dat_sub %>%
    filter(Plant_species %in% common_sp) %>%
    group_by(Plant_species) %>%
    summarise(zr = mean(zr, na.rm = TRUE), .groups = "drop") %>%
    deframe()
  
  if (is.null(tree$edge.length)) {
    message("edge.length 缺失，使用 Grafen 法赋相对分支长度")
    tree <- compute.brlen(tree, method = "Grafen")
  }
  
  phytools::phylosig(tree,
                     x      = zr_vec[tree$tip.label],
                     method = "lambda",
                     test   = TRUE)
}

cat("Overall:\n")
lambda_overall <- calc_lambda(data_zr, phylo_tree);      print(lambda_overall)

cat("\nNative + Overlap:\n")
lambda_native_ov <- calc_lambda(filter(data_zr, Plant_status %in% c("Native","Overlap")),
                                phylo_tree);             print(lambda_native_ov)

cat("\nNon‑native + Overlap:\n")
lambda_nonnative_ov <- calc_lambda(filter(data_zr, Plant_status %in% c("Non-native","Overlap")),
                                   phylo_tree);          print(lambda_nonnative_ov)

