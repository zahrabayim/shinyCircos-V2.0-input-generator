# shinycircos_input_generator.R
# Author: Zahrabayim  Abdullayeva
# Description: This script processes a normalized gene expression matrix and WGCNA module data to generate 
#              BED-format files for visualization with shinyCircos-V2.0. It includes gene coordinate retrieval 
#              via Ensembl BioMart and correlation analysis of module genes.
# Date: 2025-08-07
# Load libraries
#Install required packages if not already installed.
library(readxl) #for reading excel files
library(dplyr) #for data manipulation such as filtering and summarizing 
library(stringr) #for string operations like trimming
library(biomaRt) #for accessing Biomart  database for gene annotations 
library(forcats) #for chromosome ordering
library(tibble)  # for specific functions

# 1. Reads excel file of expression matrix 
#Make sure that the input file is placed in the "data/" directory
data <- read_xlsx("data/GSE184994_norm_log2_ratios_annotated.xlsx")
print("First input file is loaded")

# 2. Identify control and treatment columns by using Regex,value = TRUE means that it returns the actual column names. 
control_cols <- grep("^control", colnames(data), value = TRUE)
treatment_cols <- grep("^A151-", colnames(data), value = TRUE)

# 3. Clean gene names. mutate is a dplyr function that modifies columns and str_trim is a stringr function that removes whitespaces
# and filter function keeps rows where Gene rows are not empty and NA.
data <- data %>%
  dplyr::mutate(Gene = str_trim(Gene)) %>%
  dplyr::filter(Gene != "" & !is.na(Gene))

# 4. Calculate median expression for control and treatment.rowwise is a dplyr function that tells R to do the operation rowwise.
#mutate adds new columns median_control_expression is the median of values across all control columns. na.rm = TRUE ignores NA values.
#ungroup resets the data frame to the normal format.
data <- data %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    median_control_expression = median(c_across(all_of(control_cols)), na.rm = TRUE),
    median_treatment_expression = median(c_across(all_of(treatment_cols)), na.rm = TRUE)
  ) %>%
  dplyr::ungroup()

# 5. Collapse duplicates by gene. group_by is a dplyr function that treats all the rows with the same gene name as a group.
#then mean of median values are calculated. summarize only takes columns that are stated, so data_summary only has 3 columns.
data_summary <- data %>%
 dplyr:: group_by(Gene) %>%
  dplyr::summarise(
    median_control_expression = mean(median_control_expression, na.rm = TRUE),
    median_treatment_expression = mean(median_treatment_expression, na.rm = TRUE),
    .groups = "drop"
  )

# 6. Collapse full control replicate values by gene. This one calculates means of control expressions across all columns per gene.
data_controls_collapsed <- data %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(across(all_of(control_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# 7. Connect to Ensembl (mouse). Ensembl is a genome browser, and biomart is an interface that is used to connect to musculus genome.
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# 8. Retrieve coordinates
gene_locations <- getBM(
  attributes = c("external_gene_name", "chromosome_name", "start_position", "end_position"),
  filters = "external_gene_name",
  values = data_summary$Gene,
  mart = ensembl
)

# 9. Select one transcript per gene (lowest start)
gene_locations_unique <- gene_locations %>%
  dplyr::group_by(external_gene_name) %>%
  dplyr::slice_min(start_position, with_ties = FALSE) %>%
  dplyr::ungroup()

# 10. Join with expression values. joined data_summary and gene_locations_unique and added "chr" in front of chromosome numbers
joined_data <- data_summary %>%
  dplyr::inner_join(gene_locations_unique, by = c("Gene" = "external_gene_name")) %>%
  dplyr::mutate(chr = paste0("chr", chromosome_name))

# 11. Define chromosome ordering
chrom_order <- c(paste0("chr", 1:19), "chrX", "chrY", "chrMT")

# 12. Create chromosome-wise BED file (one row per chromosome, using min start and max end)
chromosome_ranges <- joined_data %>%
  dplyr::mutate(chr = factor(chr, levels = chrom_order)) %>%
  dplyr::group_by(chr) %>%
  dplyr::summarise(
    start = min(start_position, na.rm = TRUE),
    end = max(end_position, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(chr)

write.table(chromosome_ranges, "gene_coordinates.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# 13. Output with gene names
output_with_names <- joined_data %>%
  dplyr::select(chr, start = start_position, end = end_position, gene = Gene) %>%
  dplyr::mutate(chr = factor(chr, levels = chrom_order)) %>%
  dplyr::arrange(chr, start)

output_control <- joined_data %>%
  dplyr::select(chr, start = start_position, end = end_position, median_control_expression) %>%
  dplyr::mutate(chr = factor(chr, levels = chrom_order)) %>%
  dplyr::arrange(chr, start)

output_treatment <- joined_data %>%
  dplyr::select(chr, start = start_position, end = end_position, median_treatment_expression) %>%
  dplyr::mutate(chr = factor(chr, levels = chrom_order)) %>%
  dplyr::arrange(chr, start)

write.table(output_with_names, "gene_coordinates_names.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(output_control, "gene_control.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(output_treatment, "gene_treatment.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
print("Gene coordinates and expression bed files are generated")
# 14. Load WGCNA module genes
#Make sure that the input file is placed in the "data/" directory
wgcna_modules <- read_xlsx("data/1-s2.0-S2162253123000240-mmc2.xlsx", skip = 1)
print("Second input file is loaded")
red_genes <- wgcna_modules$Red %>% na.omit() %>% str_trim()
turq_genes <- wgcna_modules$Turquoise %>% na.omit() %>% str_trim()
combined_genes <- unique(c(red_genes, turq_genes))

# 15. Output coordinates + gene names for Red + Turquoise
red_turq_output <- joined_data %>%
  dplyr::filter(Gene %in% combined_genes) %>%
  dplyr::select(chr, start = start_position, end = end_position, gene = Gene) %>%
  dplyr::mutate(chr = factor(chr, levels = chrom_order)) %>%
  dplyr::arrange(chr, start)

write.table(red_turq_output, "gene_red_turquoise.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# 16. Red + Turq: Median expression values
combined_control <- joined_data %>%
  dplyr::filter(Gene %in% combined_genes) %>%
  dplyr::select(chr, start = start_position, end = end_position, value = median_control_expression) %>%
  dplyr::mutate(chr = factor(chr, levels = chrom_order)) %>%
  dplyr::arrange(chr, start)

combined_treatment <- joined_data %>%
  dplyr::filter(Gene %in% combined_genes) %>%
  dplyr::select(chr, start = start_position, end = end_position, value = median_treatment_expression) %>%
  dplyr::mutate(chr = factor(chr, levels = chrom_order)) %>%
  dplyr::arrange(chr, start)

write.table(combined_control, "gene_red_turq_control.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(combined_treatment, "gene_red_turq_treatment.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# 17. Full control replicate matrix for Red + Turq genes
rep_matrix <- data_controls_collapsed %>%
  dplyr::inner_join(gene_locations_unique, by = c("Gene" = "external_gene_name")) %>%
  dplyr::filter(Gene %in% combined_genes) %>%
  dplyr::mutate(chr = paste0("chr", chromosome_name)) %>%
  dplyr::select(chr, start = start_position, end = end_position, all_of(control_cols)) %>%
  dplyr::mutate(chr = factor(chr, levels = chrom_order)) %>%
  dplyr::arrange(chr, start)

write.table(rep_matrix, "gene_red_turq_control_replicates.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# 18. Correlation between Red + Turquoise genes (using control replicates)
expr_matrix <- data_controls_collapsed %>%
  dplyr::filter(Gene %in% combined_genes) %>%
  column_to_rownames("Gene") %>%
  dplyr::select(all_of(control_cols)) %>%
  as.matrix()

#Pearson correlation
cor_matrix <- cor(t(expr_matrix), use = "pairwise.complete.obs", method = "pearson")

cor_df <- as.data.frame(as.table(cor_matrix))
colnames(cor_df) <- c("Gene1", "Gene2", "correlation")

# Convert factors to characters before filtering
cor_df <- cor_df %>%
  dplyr::mutate(
    Gene1 = as.character(Gene1),
    Gene2 = as.character(Gene2)
  ) %>%
  dplyr::filter(Gene1 < Gene2)

coord_lookup <- gene_locations_unique %>%
  dplyr::filter(external_gene_name %in% combined_genes) %>%
  dplyr::group_by(external_gene_name) %>%
  dplyr::slice_min(start_position, with_ties = FALSE) %>%
  dplyr::mutate(chr = paste0("chr", chromosome_name)) %>%
  dplyr::select(Gene = external_gene_name, chr, start = start_position, end = end_position)

cor_with_coords <- cor_df %>%
  dplyr::left_join(coord_lookup, by = c("Gene1" = "Gene")) %>%
  dplyr::rename(chr1 = chr, start1 = start, end1 = end) %>%
  dplyr::left_join(coord_lookup, by = c("Gene2" = "Gene")) %>%
  dplyr::rename(chr2 = chr, start2 = start, end2 = end)

final_correlation_bed <- cor_with_coords %>%
  dplyr::filter(!is.na(chr1) & !is.na(chr2)) %>%
  dplyr::select(chr1, start1, end1, chr2, start2, end2, correlation)
#Writing to the file
write.table(final_correlation_bed, "gene_red_turq_correlation_pairs.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
print("Correlation file is generated")
print("script finished")