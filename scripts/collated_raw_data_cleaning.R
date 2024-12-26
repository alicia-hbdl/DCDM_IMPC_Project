library(tidyverse)

# Set working directory and load data
setwd("~/Desktop/working_directory/DCDM_project/data")
data <- read.csv("analysis_table.txt", header = TRUE, stringsAsFactors = FALSE)
sop <- read_csv("metadata/IMPC_SOP.csv", show_col_types = FALSE)

# Log file
log_file <- "collated_validation_log.txt"
file.create(log_file)

# Custom logging function
log_message <- function(message) {
  cat(message, "\n", file = log_file, append = TRUE)
  message(message) # Optionally display in the console
}

# Start logging
log_message("=== Data Validation and Transformation Log ===")

### Validate 'analysis_id'

invalid_ids <- data %>%
  filter(nchar(analysis_id) != 15 | duplicated(analysis_id))
log_message(if (nrow(invalid_ids) == 0) {
  "All analysis_id values are valid."
} else {
  paste("Invalid analysis_id values found:", nrow(invalid_ids))
})

### Validate 'gene_accession_id'

invalid_gene_lengths <- data %>% filter(nchar(gene_accession_id) < 9 | nchar(gene_accession_id) > 11)
duplicate_genes <- data %>% group_by(gene_accession_id) %>% filter(n() > 1) %>% ungroup()
log_message(if (nrow(invalid_gene_lengths) == 0) "All gene_accession_id values are valid." else paste("Invalid gene_accession_id lengths:", nrow(invalid_gene_lengths)))
log_message(if (nrow(duplicate_genes) == 0) "All gene_accession_id values are unique." else paste("Duplicate gene_accession_id values:", nrow(duplicate_genes)))

### Validate & format 'gene_symbol'

# Convert all gene_symbols to title format
data <- data %>% mutate(gene_symbol = str_to_title(gene_symbol))
# Check length of gene_symbol values
invalid_gene_symbols <- data %>% filter(nchar(gene_symbol) < 1 | nchar(gene_symbol) > 13)
log_message(if (nrow(invalid_gene_symbols) == 0) "All gene_symbol values are valid." else paste("Invalid gene_symbol values:", nrow(invalid_gene_symbols)))

### Check `gene_accession_id` and `gene_symbol` correlation

gene_symbol_correlation <- data %>%
  group_by(gene_accession_id) %>%
  filter(n_distinct(gene_symbol) > 1)
log_message(if (nrow(gene_symbol_correlation) == 0) "Perfect correlation between gene_accession_id and gene_symbol." else "Correlation issues detected.")

#### Validate and fix `mouse_strain`

valid_strains <- c("C57BL", "B6J", "C3H", "129SV")
data <- data %>%
  mutate(mouse_strain = ifelse(mouse_strain %in% valid_strains, mouse_strain, "C57BL"))
log_message("Updated mouse_strain values.") # Selected 'C57BL' after looking at unique results - there appear to be typos of C57BL

### Validate 'mouse_life_stage'

allowed_life_stages <- c("E12.5", "E15.5", "E18.5", "E9.5", "Early adult", "Late adult", "Middle aged adult")
invalid_life_stages <- data %>% filter(!mouse_life_stage %in% allowed_life_stages)
log_message(if (nrow(invalid_life_stages) == 0) "All mouse_life_stage values are valid." else paste("Invalid mouse_life_stage values:", nrow(invalid_life_stages)))

### Validate 'parameter_id'

invalid_parameter_ids <- data %>% filter(nchar(parameter_id) < 15 | nchar(parameter_id) > 18)
log_message(if (nrow(invalid_parameter_ids) == 0) {"All parameter_id values are valid."} else {paste("Invalid parameter_id values:", nrow(invalid_parameter_ids))})

### Validate 'parameter_name'

invalid_parameter_names <- data %>% filter(nchar(parameter_name) < 2 | nchar(parameter_name) > 74)
log_message(if (nrow(invalid_parameter_names) == 0) {"All parameter_name values have valid lengths."} else {paste("Invalid parameter_name values:", nrow(invalid_parameter_names))})

### Correlation between 'parameter_id' & 'parameter_name'

parameter_correlation <- data %>% group_by(parameter_id) %>% summarise(unique_names = n_distinct(parameter_name)) %>% filter(unique_names > 1)
log_message(if (nrow(parameter_correlation) == 0) {"Perfect correlation between parameter_id and parameter_name."} else {"Correlation issues detected for parameter_id and parameter_name."})

### Validate & fix 'pvalue'

data <- data %>%
  mutate(pvalue = ifelse(pvalue > 1, 1, pvalue)) # Viewing pvalues before transforming shows values >1
invalid_pvalues <- data %>% filter(pvalue < 0 | pvalue > 1 | is.na(pvalue))
log_message(if (nrow(invalid_pvalues) == 0) "All pvalue values are valid." else paste("Invalid pvalue values:", nrow(invalid_pvalues)))

### SAVE CLEANED DATA ###

write.csv(data, "cleaned_analysis_table.txt", row.names = FALSE)
log_message("Cleaned dataset saved as 'cleaned_analysis_table.txt'.")
log_message("=== End of Log ===")
