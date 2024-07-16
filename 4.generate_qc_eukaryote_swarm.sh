#!/bin/bash
# Author: Diana Proctor
# Revised July 8, 2024
# Purpose: Take the concoct bins, run eukcc, quast, and infernal on them 

# Define the number of threads based on the SLURM job allocation
THREADS=${SLURM_CPUS_ON_NODE:-1}

# Define directories
BASE_DIR="/data/proctordm/ESKAPE"
REFINEMENT_DIR="$BASE_DIR/data/03_refinement"
MEGAHIT_BINS="$REFINEMENT_DIR/03_megahit/**/concoct_bins"
SPADES_BINS="$REFINEMENT_DIR/03_spades/**/concoct_bins"
QC_EUKARYOTE_DIR="$BASE_DIR/data/04_QC/eukaryote_all"
QC_EUKARYOTE_MEGAHIT_DIR="$QC_EUKARYOTE_DIR/megahit"
QC_EUKARYOTE_SPADES_DIR="$QC_EUKARYOTE_DIR/spades"
EUKCC_DB_DIR="$BASE_DIR/data/databases/eukcc2_db_ver_1.1"
EUKCC_OUTPUT_DIR="$BASE_DIR/data/04_QC/eukaryote_all/eukcc_out"
EUKCC_INPUT_DIR="$BASE_DIR/data/04_QC/eukaryote_all/input"
EUKCC_OUTPUT_FILE="$EUKCC_OUTPUT_DIR/eukcc.csv"
EUKARYOTIC_PASSES_DIR="$BASE_DIR/data/04_QC/eukaryotic_passes"
QUAST_EUKARYOTIC_OUTPUT="$BASE_DIR/data/04_QC/eukaryote_all/quast_out"
INFERNAL_EUKARYOTIC_OUTPUT="$BASE_DIR/data/04_QC/eukaryote_all/infernal_out"
RFAM_CM_FILE="/data/proctordm/ESKAPE/data/databases/Rfam.cm"

# Ensure the necessary directories exist
directories=(
    "$BASE_DIR/data"
    "$BASE_DIR/data/04_QC"
    "$QC_EUKARYOTE_DIR"
    "$QC_EUKARYOTE_MEGAHIT_DIR"
    "$QC_EUKARYOTE_SPADES_DIR"
    "$EUKCC_DB_DIR"
    "$EUKCC_OUTPUT_DIR"
    "$EUKARYOTIC_PASSES_DIR"
    "$QUAST_EUKARYOTIC_OUTPUT"
    "$INFERNAL_EUKARYOTIC_OUTPUT"
)

for dir in "${directories[@]}"; do
    mkdir -p "$dir"
    chmod -R 755 "$dir"
done

# Check if the Rfam.cm file exists
if [ ! -f "$RFAM_CM_FILE" ]; then
    echo "Error: Rfam.cm file not found at $RFAM_CM_FILE."
    
    # URL to download Rfam.cm.gz
    URL="https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz"
    
    # Download the file
    echo "Downloading Rfam.cm.gz from $URL..."
    wget "$URL" -O "${RFAM_CM_FILE}.gz"
    
    # Check if the download was successful
    if [ ! -f "${RFAM_CM_FILE}.gz" ]; then
        echo "Download failed."
        exit 1
    fi

    # Decompress the file
    echo "Decompressing Rfam.cm.gz..."
    gunzip "${RFAM_CM_FILE}.gz"
    
    echo "File downloaded and ready to use at $RFAM_CM_FILE"
else
    echo "Rfam.cm file exists at $RFAM_CM_FILE."
fi

# Check if the EukCC database is set up
if [ ! -d "$EUKCC_DB_DIR" ]; then
    echo "EukCC database not found at $EUKCC_DB_DIR. Downloading database..."
    eukcc download_db --db_dir "$EUKCC_DB_DIR"
    if [ $? -ne 0 ]; then
        echo "Failed to download EukCC database."
        exit 1
    fi
    echo "EukCC database downloaded and set up at $EUKCC_DB_DIR."
else
    echo "EukCC database exists at $EUKCC_DB_DIR."
fi

# Function to copy and rename eukaryotic bins
copy_and_rename_eukaryotic_bins() {
    local bins_dir=$1
    local target_dir=$2
    local prefix=$3
    local copy_log="$QC_EUKARYOTE_DIR/copy_log.txt"

    if [ -d "$bins_dir" ]; then
        for bin_file in "$bins_dir"/*.fa; do
            if [ -f "$bin_file" ]; then
                bin_name=$(basename "$bin_file")
                new_bin_name="${prefix}_${bin_name}"
                if cp "$bin_file" "$target_dir/$new_bin_name"; then
                    echo "Successfully copied $bin_file to $target_dir/$new_bin_name" >> "$copy_log"
                else
                    echo "Failed to copy $bin_file to $target_dir/$new_bin_name" >> "$copy_log"
                fi
            else
                echo "Missing file: $bin_file not found." >> "$copy_log"
            fi
        done
        echo "Bins copied and renamed from $bins_dir." >> "$copy_log"
    else
        echo "No bins found in $bins_dir to copy." >> "$copy_log"
    fi
}

# Initialize the copy log file for eukaryotic bins
copy_log="$QC_EUKARYOTE_DIR/copy_log.txt"
echo "Copy log for eukaryotic bins:" > "$copy_log"

# Define a log file for debugging output
DEBUG_LOG="$REFINEMENT_DIR/debug_log.txt"

# Initialize debug log
echo "Starting the copy and rename process for eukaryotic bins" > "$DEBUG_LOG"

# Copy and rename eukaryotic bins from MegaHit
echo "Processing MegaHit refinement directories..." >> "$DEBUG_LOG"
for refined_dir in $MEGAHIT_BINS; do
    if [ -d "$refined_dir" ]; then
        echo "Processing directory: $refined_dir" >> "$DEBUG_LOG"
        sample_name=$(basename "$(dirname "$refined_dir")")
        copy_and_rename_eukaryotic_bins "$refined_dir" "$QC_EUKARYOTE_MEGAHIT_DIR" "${sample_name}_megahit_concoct_bins"
        echo "Finished processing $sample_name" >> "$DEBUG_LOG"
    else
        echo "Error: Directory does not exist: $refined_dir" >> "$DEBUG_LOG"
    fi
done

# Copy and rename eukaryotic bins from SPAdes
echo "Processing SPAdes refinement directories..." >> "$DEBUG_LOG"
for refined_dir in $SPADES_BINS; do
    if [ -d "$refined_dir" ]; then
        echo "Processing directory: $refined_dir" >> "$DEBUG_LOG"
        sample_name=$(basename "$(dirname "$refined_dir")")
        copy_and_rename_eukaryotic_bins "$refined_dir" "$QC_EUKARYOTE_SPADES_DIR" "${sample_name}_spades_concoct_bins"
        echo "Finished processing $sample_name" >> "$DEBUG_LOG"
    else
        echo "Error: Directory does not exist: $refined_dir" >> "$DEBUG_LOG"
    fi
done

echo "Copy and rename process completed." >> "$DEBUG_LOG"

# Function to run EukCC on a folder
run_eukcc() {
    local input_dir=$1
    local output_dir=$2

    mkdir -p "$output_dir"
    eukcc folder "$input_dir" --out "$output_dir" --db "$EUKCC_DB_DIR" --threads "$THREADS"
    if [ -f "$output_dir/results/eukcc.csv" ]; then
        cp "$output_dir/results/eukcc.csv" "$EUKCC_OUTPUT_FILE"
    else
        echo "EukCC results file not found in $output_dir/results" >> "$copy_log"
    fi
}

# Activate conda environments
source myconda
conda activate eukcc

echo "Running EukCC on eukaryotic bins folder..."
run_eukcc "$QC_EUKARYOTE_MEGAHIT_DIR" "$EUKCC_OUTPUT_DIR/megahit"
run_eukcc "$QC_EUKARYOTE_SPADES_DIR" "$EUKCC_OUTPUT_DIR/spades"


# Function to create QUAST Swarm files for eukaryotic bins
create_quast_swarm_eukaryote() {
    local assemblies=("$@")
    for assembly_file in "${assemblies[@]}"; do
        output_dir="$QUAST_EUKARYOTIC_OUTPUT/$(basename "${assembly_file%.*}")"
        swarm_file="$QUAST_EUKARYOTIC_OUTPUT/$(basename "${assembly_file%.*}").swarm"
        echo "metaquast.py $assembly_file -o $output_dir --max-ref-number 0" > "$swarm_file"
    done
}

# Function to create Infernal Swarm files for eukaryotic bins
create_infernal_swarm_eukaryote() {
    local bins=("$@")
    for bin_file in "${bins[@]}"; do
        swarm_file="$INFERNAL_EUKARYOTIC_OUTPUT/$(basename "$bin_file").swarm"
        echo "cmsearch -Z 1000 --hmmonly --cut_ga --noali --tblout=$INFERNAL_EUKARYOTIC_OUTPUT/$(basename "$bin_file").infernal $RFAM_CM_FILE $bin_file" > "$swarm_file"
    done
}

# Example usage: List of assemblies and bins to process
megahit_bins=("$QC_EUKARYOTE_MEGAHIT_DIR"/*.fa)
spades_bins=("$QC_EUKARYOTE_SPADES_DIR"/*.fa)

# Create QUAST Swarm files for eukaryotic bins
echo "Creating QUAST Swarm files for eukaryotic bins..."
create_quast_swarm_eukaryote "${megahit_bins[@]}"
create_quast_swarm_eukaryote "${spades_bins[@]}"

# Create Infernal Swarm files for eukaryotic bins
echo "Creating Infernal Swarm files for eukaryotic bins..."
create_infernal_swarm_eukaryote "${megahit_bins[@]}"
create_infernal_swarm_eukaryote "${spades_bins[@]}"

# Function to execute Swarm files
execute_swarm_files() {
    local swarm_dir=$1
    for swarm_file in "$swarm_dir"/*.swarm; do
        swarm -f "$swarm_file" -g 16 -t $THREADS --module quast  # Adjust module as needed
    done
}

# Execute QUAST Swarm files for eukaryotic bins
echo "Executing QUAST Swarm files for eukaryotic bins..."
execute_swarm_files "$QUAST_EUKARYOTIC_OUTPUT"

# Execute Infernal Swarm files for eukaryotic bins
echo "Executing Infernal Swarm files for eukaryotic bins..."
execute_swarm_files "$INFERNAL_EUKARYOTIC_OUTPUT"

echo "Analysis completed for eukaryotic bins."

# All operations completed successfully
echo "All operations completed successfully."
