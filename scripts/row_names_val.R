library(tidyverse)

# Set paths
setwd("~/Desktop/working_directory/DCDM_project/data")
data_dir <- "updated_data"  # Directory containing the updated data files
log_file <- "row_name_validation.log"  # Log file for results

# Load SOP and normalize row names
expected_row_names <- read_csv("metadata/IMPC_SOP.csv", show_col_types = FALSE)$dataField

# Clear the log file
file.create(log_file)

# Function to validate a file against the SOP layout and column count
validate_file <- function(file_path) {
  # Read the file
  data <- tryCatch(
    read_csv(file_path, col_names = FALSE, show_col_types = FALSE),
    error = function(e) return(NULL))
  
  # Skip files that could not be read
  if (is.null(data)) return(FALSE)
  
  # Check column count
  if (ncol(data) != 2) {
    return(list(valid = FALSE, message = "Incorrect column count"))
  }
  
  # Normalize row names and validate against SOP
  row_names <- tolower(data[[1]])
  missing_rows <- setdiff(expected_row_names, row_names)
  extra_rows <- setdiff(row_names, expected_row_names)
  
  # Create validation message
  if (length(missing_rows) > 0 || length(extra_rows) > 0) {
    message <- c(
      if (length(missing_rows) > 0) paste("Missing row names:", paste(missing_rows, collapse = ", ")),
      if (length(extra_rows) > 0) paste("Extra row names:", paste(extra_rows, collapse = ", "))
    )
    return(list(valid = FALSE, message = paste(message, collapse = "\n")))
  }
  
  list(valid = TRUE, message = "File is valid")
}

# Process all files
validation_results <- lapply(list.files(data_dir, pattern = "\\.csv$", full.names = TRUE), function(file_path) {
  result <- validate_file(file_path)
  log_message <- paste(basename(file_path), ":", result$message, "\n")
  cat(log_message, file = log_file, append = TRUE)
  result$valid
})

# Summarize and log results
valid_files <- sum(unlist(validation_results))
summary_message <- sprintf(
  "Validation Summary:\nTotal files checked: %d\nTotal valid files: %d\n", 
  length(validation_results), valid_files
)
cat(summary_message, file = log_file, append = TRUE)
message(summary_message)

