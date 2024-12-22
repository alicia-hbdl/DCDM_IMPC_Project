#!/bin/bash
#SBATCH -p msc_appbio             # Partition name
#SBATCH --ntasks=2                # Number of tasks
#SBATCH -J Rscript                # Job name

# Load the R module (if required on your HPC system)
module load r/4.3.0-gcc-13.2.0-withx-rmath-standalone-python-3.11.6

# Install R packages with a specified CRAN mirror
Rscript -e "install.packages('dplyr', repos='https://cran.r-project.org')"
Rscript -e "install.packages('readr', repos='https://cran.r-project.org')"

# Run the R script
Rscript row_names_val.R
