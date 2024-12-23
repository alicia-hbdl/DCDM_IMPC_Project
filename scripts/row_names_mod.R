library(tidyverse)

# Set paths to directories
setwd("~/Desktop/working_directory/DCDM_project/data")
data_dir <- "raw_data"
output_dir <- "updated_data"
log_file <- "row_name_modification.log"

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE)

# Load SOP file and normalize row names to lowercase
expected_row_names <- read_csv("metadata/IMPC_SOP.csv", show_col_types = FALSE)$dataField

# Clear the log file
file.create(log_file)

# Function to validate, update, and reorder row names
process_file <- function(file_path) {
  # Read the CSV file as headerless
  data <- tryCatch(
    read_csv(file_path, col_names = FALSE, show_col_types = FALSE), 
    error = function(e) return(NULL))
  
  # If file reading fails, skip processing
  if (is.null(data)) return(FALSE)
  
  # Normalize row names
  data[[1]] <- tolower(data[[1]])
  
  # Identify missing rows and add them, input NA value for second column
  missing_rows <- setdiff(expected_row_names, data[[1]])
  if (length(missing_rows) > 0) {
    missing_data <- tibble(!!colnames(data)[1] := missing_rows)
    for (i in 2:ncol(data)) {
      missing_data[[colnames(data)[i]]] <- NA
    }
    data <- bind_rows(data, missing_data)
  }
  
  # Reorder rows based on expected_row_names
  data <- arrange(data, match(data[[1]], expected_row_names))
  
  # Save the updated file as headerless
  write_csv(data, file.path(output_dir, basename(file_path)), col_names = FALSE)
 
  return(TRUE)
}

# Process files and log results
validation_results <- sapply(list.files(data_dir, pattern = "\\.csv$", full.names = TRUE), process_file)

# Log summary
summary_message <- sprintf("Validation Summary:\nTotal files checked: %d\nTotal updated files: %d\n",
                           length(validation_results), sum(validation_results))
write(summary_message, log_file)
message(summary_message)

