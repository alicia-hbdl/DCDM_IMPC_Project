library(dplyr)
library(readr)

# Set paths to directories
setwd("/scratch_tmp/grp/msc_appbio/DCDM_group3")
data_dir <- "data/raw_data"
output_dir <- "data/updated_data"
log_file <- "data/row_name_modification.log"

# Creates output directory if it doesn't already exist & supresses warnings if it does exist
dir.create(output_dir, showWarnings = FALSE)

# Load the SOP file and normalize row names to lowercase
expected_row_names <- tolower(read_csv("data/metadata/IMPC_SOP.csv", show_col_types = FALSE)$dataField)

# Clear the log file if it exists or creates it if it doesn't
file.create(log_file)

# Function to validate, update, and reorder row names
process_file <- function(file_path) {
  # Read the file specified by file_path
  data <- tryCatch(
    read_csv(file_path, col_names = TRUE, show_col_types = FALSE),
    error = function(e) {
      message(paste("Error reading file:", file_path, "-", e))
      return(NULL)
    }
  )
  if (is.null(data)) return(FALSE)
  
  # Normalize row names and identify missing rows
  data[[1]] <- tolower(data[[1]])
  missing_rows <- setdiff(expected_row_names, data[[1]])
  
  # Add missing rows with NA values in the second column
  if (length(missing_rows) > 0) {
    missing_data <- tibble(
      !!colnames(data)[1] := missing_rows,  # First column: missing row names
      !!colnames(data)[2] := NA            # Second column: NA values
    )
    data <- bind_rows(data, missing_data)
  }
  
  # Reorder rows and save the updated file
  data <- data %>% arrange(match(data[[1]], expected_row_names))
  write_csv(data, file.path(output_dir, basename(file_path)))
  
  return(TRUE)
}

# Process all files and log results
files <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)
validation_results <- sapply(files, process_file)

# Write summary to log
summary_message <- sprintf(
  "Validation Summary:\nTotal files checked: %d\nTotal updated files: %d\n",
  length(validation_results), sum(validation_results)
)
cat(summary_message, file = log_file)
message(summary_message)

