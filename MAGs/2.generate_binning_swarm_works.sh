#!/bin/bash

# Error handling
error_exit() {
    echo "$1" 1>&2
    exit 1
}

# Define paths to input
DATA_DIR="/data/proctordm/ESKAPE/data"
BASE_DIR="$DATA_DIR"
READS_DIR="/data/proctordm/ESKAPE/00_reads"

BINNING_DIR="$DATA_DIR/02_binning"
SCRIPTS_DIR="$BINNING_DIR/scripts"
LOG_DIR="$BINNING_DIR/logs"
CONTAINER_DIR="/data/proctordm/ESKAPE/container"
SIF_FILE="$CONTAINER_DIR/MAG_wf_containers_metawrap.sif"

MEGAHIT_CONTIGS_LIST="$BASE_DIR/01_assembly/01_assembly_megahit/final_contigs_list.txt"
SPADES_CONTIGS_LIST="$BASE_DIR/01_assembly/01_assembly_spades/final_contigs_list.txt"
MIN_LENGTH=5000
# Define resources
THREADS=32
MEMORY=128
TIME="24:00:00"

BINNING_CMD="singularity run $SIF_FILE metawrap binning --metabat2 --maxbin2 --concoct -l $MIN_LENGTH -t $THREADS -m $MEMORY -a \${ASSEMBLY} -o \${OUTPUT} \${READ1} \${READ2}"

# Ensure directories exist
mkdir -p "$SCRIPTS_DIR"
mkdir -p "$LOG_DIR"
mkdir -p "$BINNING_DIR/megahit"
mkdir -p "$BINNING_DIR/spades"

# Create the binning.sh script
cat << EOF > "$SCRIPTS_DIR/binning.sh"
#!/usr/bin/bash
ASSEMBLY=\$1
OUTPUT=\$2
READ1=\$3 
READ2=\$4
module load singularity
source /usr/local/current/singularity/app_conf/sing_binds
$BINNING_CMD
EOF

# Make the script executable
chmod +x "$SCRIPTS_DIR/binning.sh"

# Read assemblies into arrays
mapfile -t MEGAHIT_ASSEMBLIES < "$MEGAHIT_CONTIGS_LIST"
mapfile -t SPADES_ASSEMBLIES < "$SPADES_CONTIGS_LIST"

# Create lists of forward and reverse reads
ls -d $READS_DIR/*_1.fastq > READ1.list
ls -d $READS_DIR/*_2.fastq > READ2.list

# Check if READ1.list and READ2.list were created successfully
if [[ -f "READ1.list" && -f "READ2.list" ]]; then
    echo "READ1.list and READ2.list created successfully."
else
    error_exit "Failed to create READ1.list and/or READ2.list."
fi

# Read lists into arrays
echo "Reading lists into arrays..."
mapfile -t READ1 < READ1.list
mapfile -t READ2 < READ2.list

# Extract sample names from reads
read_samples=()
for read in "${READ1[@]}"; do
    read_samples+=($(basename "$read" | sed 's/_1.fastq//'))
done

# Extract sample names from assemblies
megahit_samples=()
for assembly in "${MEGAHIT_ASSEMBLIES[@]}"; do
    megahit_samples+=($(basename "$(dirname "$assembly")" | sed 's/_megahit_out//'))
done

spades_samples=()
for assembly in "${SPADES_ASSEMBLIES[@]}"; do
    spades_samples+=($(basename "$(dirname "$assembly")" | sed 's/_spades_out//'))
done

# Track mismatches separately for MegaHit and Spades
megahit_matching_samples=()
spades_matching_samples=()

for sample in "${read_samples[@]}"; do
    if [[ " ${megahit_samples[@]} " =~ " ${sample} " ]]; then
        megahit_matching_samples+=("$sample")
    fi
    
    if [[ " ${spades_samples[@]} " =~ " ${sample} " ]]; then
        spades_matching_samples+=("$sample")
    fi
done

# Generate the swarm file for MegaHit matching samples
MEGAHIT_SWARM_FILE="$BINNING_DIR/megahit_swarm_file"
> "$MEGAHIT_SWARM_FILE"
for sample in "${megahit_matching_samples[@]}"; do
    read1="$READS_DIR/${sample}_1.fastq"
    read2="$READS_DIR/${sample}_2.fastq"
    megahit_assembly=$(grep "/${sample}_megahit_out" "$MEGAHIT_CONTIGS_LIST")

    if [ -f "$megahit_assembly" ] && [ -f "$read1" ] && [ -f "$read2" ]; then
        megahit_output_dir="$BINNING_DIR/megahit/${sample}_megahit_binning_out"
        mkdir -p "$megahit_output_dir"
        echo "bash $SCRIPTS_DIR/binning.sh $megahit_assembly $megahit_output_dir $read1 $read2" >> "$MEGAHIT_SWARM_FILE"
    else
        echo "Warning: Required file(s) missing or not readable for sample $sample. Skipping..."
    fi
done

# Generate the swarm file for Spades matching samples
SPADES_SWARM_FILE="$BINNING_DIR/spades_swarm_file"
> "$SPADES_SWARM_FILE"
for sample in "${spades_matching_samples[@]}"; do
    read1="$READS_DIR/${sample}_1.fastq"
    read2="$READS_DIR/${sample}_2.fastq"
    spades_assembly=$(grep "/${sample}_spades_out" "$SPADES_CONTIGS_LIST")

    if [ -f "$spades_assembly" ] && [ -f "$read1" ] && [ -f "$read2" ]; then
        spades_output_dir="$BINNING_DIR/spades/${sample}_spades_binning_out"
        mkdir -p "$spades_output_dir"
        echo "bash $SCRIPTS_DIR/binning.sh $spades_assembly $spades_output_dir $read1 $read2" >> "$SPADES_SWARM_FILE"
    else
        echo "Warning: Required file(s) missing or not readable for sample $sample. Skipping..."
    fi
done

echo "Swarm files created successfully."
echo "MegaHit swarm file: $MEGAHIT_SWARM_FILE"
echo "Spades swarm file: $SPADES_SWARM_FILE"

# Print the number of samples to be processed
echo "Number of MegaHit matching samples to be processed: ${#megahit_matching_samples[@]}"
echo "Number of Spades matching samples to be processed: ${#spades_matching_samples[@]}"

# Launch swarm jobs and print job IDs
echo "Launching MegaHit swarm job..."
megahit_jobid=$(swarm -f "$MEGAHIT_SWARM_FILE" -t $THREADS -g $MEMORY --time $TIME --logdir "$LOG_DIR")
echo "MegaHit swarm job launched with job ID: $megahit_jobid"

echo "Launching Spades swarm job..."
spades_jobid=$(swarm -f "$SPADES_SWARM_FILE" -t $THREADS -g $MEMORY --time $TIME --logdir "$LOG_DIR")
echo "Spades swarm job launched with job ID: $spades_jobid"

echo "Jobs launched successfully."

# Wait for swarm jobs to finish
echo "Waiting for MegaHit swarm jobs to finish..."
while squeue -u $USER -j $megahit_jobid > /dev/null 2>&1; do
    sleep 60
done

echo "Waiting for Spades swarm jobs to finish..."
while squeue -u $USER -j $spades_jobid > /dev/null 2>&1; do
    sleep 60
done

# Write message indicating that swarm jobs finished successfully
echo "Swarm jobs completed successfully" >> "$LOG_DIR/swarm_completion.log"

# Analyze MegaHit log files for failures and create a list of failed samples
echo "Analyzing MegaHit log files for failures..."
MEGAHIT_FAILED_SAMPLES_LIST="$BINNING_DIR/megahit_failed_samples_list"
> $MEGAHIT_FAILED_SAMPLES_LIST
for logfile in $LOG_DIR/*.o*; do
    if ! grep -q "PIPELINE SUCCESSFULLY FINISHED" "$logfile"; then
        sample=$(grep -oP '(?<=final.contigs.fa ).*(?= )' "$logfile" | head -n 1)
        if [ -n "$sample" ]; then
            echo "$sample" >> $MEGAHIT_FAILED_SAMPLES_LIST
        fi
    fi
done

# Check if MegaHit failed samples list was created successfully and write "none" if empty
if [[ -s "$MEGAHIT_FAILED_SAMPLES_LIST" ]]; then
    echo "MegaHit failed samples list created at $MEGAHIT_FAILED_SAMPLES_LIST"
else
    echo "none" > $MEGAHIT_FAILED_SAMPLES_LIST
fi

# Analyze Spades log files for failures and create a list of failed samples
echo "Analyzing Spades log files for failures..."
SPADES_FAILED_SAMPLES_LIST="$BINNING_DIR/spades_failed_samples_list"
> $SPADES_FAILED_SAMPLES_LIST
for logfile in $LOG_DIR/*.o*; do
    if ! grep -q "PIPELINE SUCCESSFULLY FINISHED" "$logfile"; then
        sample=$(grep -oP '(?<=final.contigs.fa ).*(?= )' "$logfile" | head -n 1)
        if [ -n "$sample" ]; then
            echo "$sample" >> $SPADES_FAILED_SAMPLES_LIST
        fi
    fi
done

# Check if Spades failed samples list was created successfully and write "none" if empty
if [[ -s "$SPADES_FAILED_SAMPLES_LIST" ]]; then
    echo "Spades failed samples list created at $SPADES_FAILED_SAMPLES_LIST"
else
    echo "none" > $SPADES_FAILED_SAMPLES_LIST
fi

echo "Script completed successfully."

# Function to count the number of failed samples
count_failed_samples() {
    failed_samples_file="$1"
    grep -cv '^none$' "$failed_samples_file"
}

# Get the count of failed samples for MegaHit and SPAdes
megahit_failed_count=$(count_failed_samples "$MEGAHIT_FAILED_SAMPLES_LIST")
spades_failed_count=$(count_failed_samples "$SPADES_FAILED_SAMPLES_LIST")

# Print the results
echo "Number of samples that failed MegaHit binning: $megahit_failed_count"
echo "Number of samples that failed SPAdes binning: $spades_failed_count"
