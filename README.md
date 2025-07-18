# Description of repository

The documents in this repository provide the R scripts, data, and supplementary files used for analyzing global latitudinal patterns of herbivory in native versus non-native plants, including the roles of climatic gradients, herbivore feeding guilds, and plant phylogeny. For comprehensive results and methodological details, please refer to our study titled:

> **"Herbivory increases towards lower latitudes in native but not introduced plants."**

**Corresponding author:** Rui-Ting Ju ([jurt@fudan.edu.cn](mailto:jurt@fudan.edu.cn))   

The webpage containing the neatly formatted supplementary material (code, figures, and data) can be found here:  
[https://yaolinguo.github.io/latitudinal_herbivory_gradients/](https://yaolinguo.github.io/latitudinal_herbivory_gradients/)

**Data folder:**

- `Data_zr_herbivory-latitude.csv`   This file contains Fisher’s z‑transformed effect sizes (`zr`) representing the strength and direction of latitudinal herbivory gradients, based on correlations between herbivory intensity and latitude across studies. It was used to assess whether the strength and direction of latitudinal herbivory gradients differ between native and non-native plant species.
  
- `Data_zr_herbivory-bioclimate.csv`  This dataset contains Fisher’s z‑transformed effect sizes (`zr`), calculated from correlations between herbivory intensity and 19 bioclimatic variables (Bio1–Bio19). It was used to evaluate how herbivory patterns respond to climatic gradients across latitude.

- `README.md`  describes the meaning of every column in both datasets.

**Code folder:**

- `R code 1 - Geospatial mapping and overall effect size synthesis - 0713.R`  Cleans the raw dataset, plots the geographic distribution of all study sites on a world map, estimates the overall strength of the latitudinal herbivory gradient, and contrasts that gradient between native and non‑native plant species.

- `R code 2 - Herbivory feeding guild effects - 0713.R`  Conducts subgroup meta‑analyses by herbivore feeding guilds, produces forest plots, and tests whether feeding guild differentially moderates herbivory in native versus non‑native plants.

- `R code 3 - Bioclimatic variable effects - 0713.R`  Quantifies how bioclimatic variables influence the latitudinal pattern of herbivory.

- `R code 4 - Publication bias assessment - 0713.R`  Detects and corrects publication bias using funnel plots, Egger’s regression, and trim‑and‑fill procedures, and exports comprehensive bias‑diagnostic outputs.

- `R code 5 - Phylogenetic tree construction and analysis - 0713.R`  Builds a phylogenetic tree and calculates Pagel’s λ to evaluate phylogenetic signal in the strength of the latitudinal herbivory gradient.

- `R code 6 - Latitudinal range effects - 0713.R`  Uses regression to test how a species’ latitudinal range width (geographic breadth) modulates the strength of the latitudinal herbivory gradient.

**Figures folder:**

- `Figure 1 - 0713.pdf`  Conceptual models illustrating six potential scenarios of latitudinal gradients in herbivory on native and non-native plants.

- `Figure 2 - 0713.pdf`  Map of all herbivory observation sites included in this meta-analysis, and the average strengths of latitudinal herbivory gradients.

- `Figure 3 - 0713.pdf`  The average strengths of latitudinal herbivory gradients for different herbivore feeding guilds across all plants, and separately for native and non-native plants.

- `Figure 4 - 0713.pdf`  Results from REML meta‑analyses showing the effects of bioclimatic variables on herbivory.
  
- `Figure 5 - 0713.pdf`  Tests for publication bias in the dataset.
       
- `Figure S1 - 0713.pdf`  Phylogenetic tree of the collected plant species and their average strengths of latitudinal herbivory gradients.
         
- `Figure S2 - 0713.pdf`  The relationships between herbivory effect sizes and latitudinal extension.
