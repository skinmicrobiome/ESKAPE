#!/bin/bash

# Fix eukaryotic bin names

# Directory where multifasta files are located
input_dir="/data/proctordm/ESKAPE/all_mags"

# Directory where modified fasta files will be stored
output_dir="/data/proctordm/ESKAPE/all_mags/header"
mkdir -p "$output_dir"

# Function to rename multifasta headers
rename_multifasta_headers() {
    local fasta_file="$1"
    local output_file="$2"

    # Read fasta file line by line
    while IFS= read -r line; do
        if [[ $line == ">"* ]]; then
            # Extract the initial header without ">"
            initial_header=$(echo "$line" | cut -d'>' -f2)
            
            # Extract the bin name
            bin_name=$(basename "$fasta_file" .fa)

            # Reconstruct the header with desired format
            new_header=">$bin_name $initial_header"

            echo "$new_header"
        else
            echo "$line"
        fi
    done < "$fasta_file" > "$output_file"
}

# Loop through each multifasta file in the input directory
for fasta_file in "$input_dir"/*.fa; do
    if [ -f "$fasta_file" ]; then
        # Define output file path
        output_file="$output_dir/$(basename "$fasta_file")"

        # Perform header renaming
        rename_multifasta_headers "$fasta_file" "$output_file"
    fi
done

echo "Header renaming completed."

