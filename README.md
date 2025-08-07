shinyCircos-V2.0-input-generator

##Description
This R script processes a normalized gene expression matrix along with WGCNA module gene lists (Red & Turquoise) to generate multiple BED-format files containing gene annotations, correlation values, and expression data. These files are designed for comprehensive visualization with shinyCircos-V2.0. The script also retrieves gene genomic coordinates via Ensembl BioMart and computes pairwise gene correlations within and between modules.

##Data source
The input files used in this analysis were obtained from:

Yazar, V., Yilmaz, I. C., Bulbul, A., Klinman, D. M., & Gursel, I. (2023). Gene network landscape of mouse splenocytes reveals integrin complex as the A151 ODN-responsive hub molecule in the immune transcriptome. Molecular therapy. Nucleic acids, 31, 553–565. https://doi.org/10.1016/j.omtn.2023.02.004

##How to Use

1. Clone the repository:
   git clone https://github.com/zahrabayim/shinyCircos-V2.0-input-generator.git

2. Open the R script (shinycircos_input_generator.R) in RStudio or your preferred R environment.

3. Place the input files (GSE184994_norm_log2_ratios_annotated.xlsx) and (1-s2.0-S2162253123000240-mmc2.xlsx) inside a folder named data/. GSE184994_norm_log2_ratios_annotated.xlsx is a normalized expression matrix in Excel format,and 1-s2.0-S2162253123000240-mmc2.xlsx is an Excel file that consists of gene lists for each module. 

4. Run the script to generate BED files for shinyCircos-V2.0 visualization.

##Requirements

- R version 4.0.0 or higher
- R packages: readxl, dplyr, stringr, biomaRt, forcats, tibble

To install required packages, run these commands in R:

install.packages(c("readxl", "dplyr", "stringr", "forcats", "tibble"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("biomaRt")

##Author

Zahrabayim Abdullayeva


