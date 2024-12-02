#!/bin/bash

# List of files
files=("Disease_information.txt" "IMPC_parameter_description.txt" "IMPC_procedure.txt")

# Directory containing the files
dir="/scratch_tmp/grp/msc_appbio/DCDM_group3/metadata"

# Loop through each file
for file in "${files[@]}"; do
    echo "===== ${file} ====="
    head -n 5 "${dir}/${file}" | column -t -s ','
    echo
done
