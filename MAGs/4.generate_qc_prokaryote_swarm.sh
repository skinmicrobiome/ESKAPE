#!/bin/bash
#Diana Proctor
# Purpose: Take the refined prokaryotic bins, run checkm2 and gunc on them

# Define the number of threads based on the SLURM job allocation
THREADS=${SLURM_CPUS_ON_NODE:-1}


# Define directories
BASE_DIR="/data/proctordm/ESKAPE"
REFINEMENT_DIR="$BASE_DIR/data/03_refinement"
MEGAHIT_REFINEMENT_DIR="$REFINEMENT_DIR/03_megahit"
SPADES_REFINEMENT_DIR="$REFINEMENT_DIR/03_spades"
QC_PROKARYOTE_DIR="$BASE_DIR/data/04_QC/QC_prokaryote_all"
MEGAHIT_QC_DIR="$QC_PROKARYOTE_DIR/qc_megahit"
SPADES_QC_DIR="$QC_PROKARYOTE_DIR/qc_spades"
RFAM_CM_FILE="$BASE_DIR/data/databases/Rfam.cm"
CHECKM_DB="/data/proctordm/ESKAPE/data/databases/CheckM2_database/uniref100.KO.1.dmnd"

# Ensure the necessary directories exist
directories=(
    "$MEGAHIT_QC_DIR/input"
    "$SPADES_QC_DIR/input"
    "$MEGAHIT_QC_DIR/gunc_out/logs"
    "$SPADES_QC_DIR/gunc_out/logs"
    "$MEGAHIT_QC_DIR/quast_output"
    "$SPADES_QC_DIR/quast_output"
    "$MEGAHIT_QC_DIR/infernal_output"
    "$SPADES_QC_DIR/infernal_output"
    "$BASE_DIR/data/databases"
)

for dir in "${directories[@]}"; do
    mkdir -p "$dir"
    chmod -R 755 "$dir"
done

# Download and set up CheckM2 if not already installed
if [ ! -d "$BASE_DIR/checkm2" ]; then
    git clone --recursive https://github.com/chklovski/checkm2.git "$BASE_DIR/checkm2"
    cd "$BASE_DIR/checkm2"
    conda env create -n checkm2 -f checkm2.yml
    conda activate checkm2
    conda install python=3.8
    bin/checkm2 -h
    cd "$BASE_DIR"
else
    echo "CheckM2 already set up."
fi

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

# Function to copy and rename prokaryotic bins
copy_and_rename_prokaryotic_bins() {
    local bins_dir=$1
    local target_dir=$2
    local prefix=$3

    echo "Copying prokaryotic bins from $bins_dir to $target_dir with prefix $prefix..."
    if [ -d "$bins_dir" ]; then
        for bin_file in "$bins_dir"/*.fa; do
            bin_name=$(basename "$bin_file")
            new_bin_name="${prefix}_${bin_name}"
            echo "Copying $bin_file to $target_dir/$new_bin_name"
            cp "$bin_file" "$target_dir/$new_bin_name"
        done
        echo "Prokaryotic bins copied and renamed from $bins_dir." >> "$copy_log"
    else
        echo "No prokaryotic bins found in $bins_dir to copy." >> "$copy_log"
    fi
}

# Function to create a Slurm job script for running GUNC
run_gunc() {
    local input_dir=$1
    local output_dir=$2
    local gunc_job_file="$output_dir/run_gunc.sh"

    echo "Creating GUNC Slurm job script at $gunc_job_file..."
    cat << EOF > "$gunc_job_file"
#!/bin/bash
#SBATCH --job-name=gunc_analysis
#SBATCH --output=$output_dir/logs/gunc_%j.out
#SBATCH --error=$output_dir/logs/gunc_%j.err
#SBATCH --cpus-per-task=$THREADS
#SBATCH --mem=128G
#SBATCH --time=24:00:00

source myconda
conda activate gunc

gunc run --input_dir "$input_dir" -r "$GUNC_DB_DIR" --out_dir "$output_dir" --threads $THREADS
EOF

    chmod +x "$gunc_job_file"
    sbatch "$gunc_job_file"
}

# Function to create CheckM2 Swarm files
create_checkm2_swarm() {
    local input_dir=$1
    local output_dir=$2
    local swarm_file=$3

    echo "Creating CheckM2 Swarm file for $input_dir..."
    mkdir -p "$output_dir/logs"
    echo "checkm2 predict --threads $THREADS --input $input_dir --output-directory $output_dir -x fa --force --database_path $CHECKM_DB" > "$swarm_file"
}

# Function to create QUAST Swarm files
create_quast_swarm() {
    local input_dir=$1
    local output_dir=$2
    local swarm_file=$3

    echo "Creating QUAST Swarm file for $input_dir..."
    mkdir -p "$output_dir/logs"
    for assembly in "$input_dir"/*.fa; do
        echo "metaquast.py $assembly -o $output_dir/$(basename "$assembly" .fa)_quast" >> "$swarm_file"
    done
}

# Function to create Infernal Swarm files
create_infernal_swarm() {
    local input_dir=$1
    local output_dir=$2
    local swarm_file=$3

    echo "Creating Infernal Swarm file for $input_dir..."
    mkdir -p "$output_dir/logs"
    for bin in "$input_dir"/*.fa; do
        echo "cmsearch -Z 1000 --hmmonly --cut_ga --noali --tblout=$output_dir/$(basename "$bin" .fa).infernal $RFAM_CM_FILE $bin" >> "$swarm_file"
    done
}

# Initialize the copy log file for prokaryotic bins
copy_log="$QC_PROKARYOTE_DIR/copy_log.txt"
echo "Copy log for prokaryotic bins:" > "$copy_log"

# Copy and rename prokaryotic bins for Megahit
echo "Copying and renaming prokaryotic bins for Megahit..."
for metawrap_50_10_bin_dir in "$MEGAHIT_REFINEMENT_DIR"/*/metawrap_50_10_bins; do
    sample_name=$(basename "$(dirname "$metawrap_50_10_bin_dir")")
    copy_and_rename_prokaryotic_bins "$metawrap_50_10_bin_dir" "$MEGAHIT_QC_DIR/input" "${sample_name}_megahit"
done

# Copy and rename prokaryotic bins for Spades
echo "Copying and renaming prokaryotic bins for Spades..."
for metawrap_50_10_bin_dir in "$SPADES_REFINEMENT_DIR"/*/metawrap_50_10_bins; do
    sample_name=$(basename "$(dirname "$metawrap_50_10_bin_dir")")
    copy_and_rename_prokaryotic_bins "$metawrap_50_10_bin_dir" "$SPADES_QC_DIR/input" "${sample_name}_spades"
done

echo "Creating GUNC, CheckM2, QUAST, and Infernal Swarm files for Megahit and Spades..."

# Create and submit GUNC Slurm jobs for Megahit and Spades
# Activate conda environments
source myconda
conda activate gunc

run_gunc "$MEGAHIT_QC_DIR/input" "$MEGAHIT_QC_DIR/gunc_out"
run_gunc "$SPADES_QC_DIR/input" "$SPADES_QC_DIR/gunc_out"

# Create CheckM2 Swarm files for Megahit and Spades
create_checkm2_swarm "$MEGAHIT_QC_DIR/input" "$MEGAHIT_QC_DIR/checkm2_out" "$MEGAHIT_QC_DIR/checkm2.swarm"
create_checkm2_swarm "$SPADES_QC_DIR/input" "$SPADES_QC_DIR/checkm2_out" "$SPADES_QC_DIR/checkm2.swarm"

# Create QUAST Swarm files for Megahit and Spades
create_quast_swarm "$MEGAHIT_QC_DIR/input" "$MEGAHIT_QC_DIR/quast_output" "$MEGAHIT_QC_DIR/quast.swarm"
create_quast_swarm "$SPADES_QC_DIR/input" "$SPADES_QC_DIR/quast_output" "$SPADES_QC_DIR/quast.swarm"

# Create Infernal Swarm files for Megahit and Spades
create_infernal_swarm "$MEGAHIT_QC_DIR/input" "$MEGAHIT_QC_DIR/infernal_output" "$MEGAHIT_QC_DIR/infernal.swarm"
create_infernal_swarm "$SPADES_QC_DIR/input" "$SPADES_QC_DIR/infernal_output" "$SPADES_QC_DIR/infernal.swarm"

# Submit CheckM2 Swarm jobs
echo "Submitting CheckM2 Swarm jobs for Megahit and Spades..."
conda activate checkm2
swarm -f "$MEGAHIT_QC_DIR/checkm2.swarm" --job-name checkm2_megahit -t $THREADS -g 128 --time 24:00:00 --logdir "$MEGAHIT_QC_DIR/checkm2_out/logs"
swarm -f "$SPADES_QC_DIR/checkm2.swarm" --job-name checkm2_spades -t $THREADS -g 128 --time 24:00:00 --logdir "$SPADES_QC_DIR/checkm2_out/logs"

# Submit QUAST Swarm jobs
echo "Submitting QUAST Swarm jobs for Megahit and Spades..."
module load quast
swarm -f "$MEGAHIT_QC_DIR/quast.swarm" --job-name quast_megahit -t $THREADS -g 128 --time 24:00:00 --logdir "$MEGAHIT_QC_DIR/quast_output/logs"
swarm -f "$SPADES_QC_DIR/quast.swarm" --job-name quast_spades -t $THREADS -g 128 --time 24:00:00 --logdir "$SPADES_QC_DIR/quast_output/logs"

# Submit Infernal Swarm jobs
conda activate infernal
echo "Submitting Infernal Swarm jobs for Megahit and Spades..."
swarm -f "$MEGAHIT_QC_DIR/infernal.swarm" --job-name infernal_megahit -t $THREADS -g 128 --time 24:00:00 --logdir "$MEGAHIT_QC_DIR/infernal_output/logs"
swarm -f "$SPADES_QC_DIR/infernal.swarm" --job-name infernal_spades -t $THREADS -g 128 --time 24:00:00 --logdir "$SPADES_QC_DIR/infernal_output/logs"

echo "All tasks submitted. Check the Slurm job logs for details."
