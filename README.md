# Description of repository

The documents in this repository provide the R scripts, data, and supplementary files used for analyzing global latitudinal patterns of herbivory in native versus non-native plants, including the roles of climatic gradients, herbivore feeding guilds, and plant phylogeny. For comprehensive results and methodological details, please refer to our study titled:

> **"Herbivory increases towards lower latitudes in native but not introduced plants."**

**Corresponding author:** Rui-Ting Ju ([jurt@fudan.edu.cn](mailto:jurt@fudan.edu.cn))   

The webpage containing the neatly formatted supplementary material (code, figures, and data) can be found here:  
[https://yaolinguo.github.io/latitudinal_gradient_in_herbivory/](https://yaolinguo.github.io/latitudinal_gradient_in_herbivory/)

**Data folder:**

- `meta_analysis_dataset.csv`  
  The source data used for all meta-analyses

- `Eligibility_Workbook.csv`  
  Table documenting the reasons for inclusion or exclusion of each study (as presented in the Supplementary Material)

- `Outcome.descriptions.csv`  
  Table documenting how herbivory was measured across different studies (also presented in the Supplementary Material)

**Figures folder:**

- `ForestPlot_large.pdf` / `ForestPlot_large.png`  
  High-resolution version of the forest plot featured in the HTML supplementary material

**Code folder:**

- `.R functions`  
  Helpful R functions sourced outside the main analysis RMarkdown (used in `meta_analysis_herbivory_latitude.rmd`)

- `meta_analysis_herbivory_latitude.rmd`  
  The complete RMarkdown with all analysis code used in the manuscript.  
  The HTML-rendered version of this file is named `index.html` and hosted here:  
  [https://yaolinguo.github.io/latitudinal_gradient_in_herbivory/](https://yaolinguo.github.io/latitudinal_gradient_in_herbivory/)

- `Response_to_Reviewers.rmd` / `.pdf`  
  The response to peer reviewers, along with a short bibliography (`bibliography_meta_analysis.bib`) used to format the document
