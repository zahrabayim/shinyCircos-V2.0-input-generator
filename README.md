shinyCircos-V2.0-input-generator

This R script processes a normalized gene expression matrix and WGCNA module gene lists (Red & Turquoise) to generate BED-format files for visualization with shinyCircos-V2.0. It retrieves gene coordinates via Ensembl BioMart and performs correlation analysis on module genes.

How to Use

1. Clone the repository:
   git clone https://github.com/zahrabayim/shinyCircos-V2.0-input-generator.git

2. Open the R script (shinycircos_input_generator.R) in RStudio or your preferred R environment.

3. Place the input files (GSE184994_norm_log2_ratios_annotated.xlsx) and (1-s2.0-S2162253123000240-mmc2.xlsx) inside a folder named data/.

4. Run the script to generate BED files and correlation data for shinyCircos visualization.

Requirements

- R version 4.0.0 or higher
- R packages: readxl, dplyr, stringr, biomaRt, forcats, tibble

To install required packages, run these commands in R:

install.packages(c("readxl", "dplyr", "stringr", "forcats", "tibble"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("biomaRt")

Author

Zahrabayim Abdullayeva


