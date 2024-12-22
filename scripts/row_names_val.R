library(dplyr)
library(readr)

# Set paths
setwd("/scratch_tmp/grp/msc_appbio/DCDM_group3")
data_dir <- "data/updated_data"  # Directory containing the updated data files
log_file <- "data/row_name_validation.log"  # Log file for results

# Load SOP and normalize row names
expected_row_names <- tolower(read_csv("data/metadata/IMPC_SOP.csv", show_col_types = FALSE)$dataField)

# Clear the log file
file.create(log_file)

# Function to validate a file against the SOP layout and column count
validate_file <- function(file_path, expected_row_names) {
  # Read the file
  data <- tryCatch(
    read_csv(file_path, col_names = TRUE, show_col_types = FALSE),
    error = function(e) {
      message(paste("Error reading file:", file_path, "-", e))
      return(list(valid = FALSE, message = "Error reading file"))
    }
  )
  
  # Skip files that could not be read
  if (is.null(data)) return(list(valid = FALSE, message = "File could not be read"))
  
  # Check column count
  if (ncol(data) != 2) {
    return(list(valid = FALSE, message = paste("Invalid column count:", ncol(data), "(expected: 2)")))
  }
  
  # Normalize row names (first column) and check against SOP
  row_names <- tolower(data[[1]])
  missing_rows <- setdiff(expected_row_names, row_names)
  extra_rows <- setdiff(row_names, expected_row_names)
  
  # Validate row names
  if (length(missing_rows) > 0 || length(extra_rows) > 0) {
    message <- ""
    if (length(missing_rows) > 0) {
      message <- paste(message, "Missing row names:", paste(missing_rows, collapse = ", "), "\n")
    }
    if (length(extra_rows) > 0) {
      message <- paste(message, "Extra row names:", paste(extra_rows, collapse = ", "), "\n")
    }
    return(list(valid = FALSE, message = message))
  }
  
  return(list(valid = TRUE, message = "File is valid"))
}

# Process all files in the data directory
files <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)
validation_results <- lapply(files, function(file_path) {
  result <- validate_file(file_path, expected_row_names)
  message <- paste(basename(file_path), ":", result$message, "\n")
  cat(message, file = log_file, append = TRUE)
  return(result$valid)
})

# Summarize results
valid_files <- sum(unlist(validation_results))
total_files <- length(validation_results)

# Write summary to log
summary_message <- paste(
  "Validation Summary:\n",
  "Total files checked:", total_files, "\n",
  "Total valid files:", valid_files, "\n"
)
cat(summary_message, file = log_file, append = TRUE)
message(summary_message)
