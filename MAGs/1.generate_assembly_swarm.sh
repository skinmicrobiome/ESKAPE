#!/bin/bash
# Author: Diana Proctor
# Revised July 5, 2024

# Define paths and directories
BASE_DIR="/data/proctordm/ESKAPE"
READS_DIR="$BASE_DIR/00_reads"
MEGAHIT_ASSEMBLY_DIR="$BASE_DIR/data/01_assembly/01_assembly_megahit"
SPADES_ASSEMBLY_DIR="$BASE_DIR/data/01_assembly/01_assembly_spades"
MEGAHIT_SWARM_FILE="$MEGAHIT_ASSEMBLY_DIR/megahit.swarm"
SPADES_SWARM_FILE="$SPADES_ASSEMBLY_DIR/spades.swarm"
MEGAHIT_LOG_DIR="$BASE_DIR/data/01_assembly/megahit_swarm_logs"
SPADES_LOG_DIR="$BASE_DIR/data/01_assembly/spades_swarm_logs"

# Lists for MegaHit
MEGAHIT_FINAL_CONTIGS_LIST="$MEGAHIT_ASSEMBLY_DIR/final_contigs_list.txt"
MEGAHIT_FAILED_SAMPLES_LIST="$MEGAHIT_ASSEMBLY_DIR/failed_samples_list.txt"

# Lists for SPAdes
SPADES_FINAL_CONTIGS_LIST="$SPADES_ASSEMBLY_DIR/final_contigs_list.txt"
SPADES_FAILED_SAMPLES_LIST="$SPADES_ASSEMBLY_DIR/failed_samples_list.txt"

# Default swarm settings
THREADS=16
MEMORY=98
TIME="24:00:00"
JOB_NAME="assembly"

# Print starting message
echo "Starting the script..."

# Create directories if they do not exist
echo "Creating directories..."
mkdir -p "$MEGAHIT_ASSEMBLY_DIR"
mkdir -p "$SPADES_ASSEMBLY_DIR"
mkdir -p "$MEGAHIT_LOG_DIR"
mkdir -p "$SPADES_LOG_DIR"

# Create megahit.sh script
echo "Creating megahit.sh script..."
cat << EOF > "$MEGAHIT_ASSEMBLY_DIR/megahit.sh"
#!/usr/bin/bash

# Check if the correct number of arguments is provided
if [ "\$#" -ne 3 ]; then
    echo "Usage: \$0 <read1> <read2> <output>"
    exit 1
fi

READ1=\$1
READ2=\$2
OUTPUT=\$3

# Load the megahit module
module load megahit
if [ \$? -ne 0 ]; then
    echo "Failed to load megahit module"
    exit 1
fi

# Run megahit
megahit -1 \${READ1} -2 \${READ2} -o \${OUTPUT} --memory $MEMORY -t $THREADS
if [ \$? -ne 0 ]; then
    echo "Megahit failed for \$OUTPUT"
    exit 1
fi

echo "Megahit completed successfully for \$OUTPUT"
EOF
# Check if megahit.sh was created successfully
if [[ -f "$MEGAHIT_ASSEMBLY_DIR/megahit.sh" ]]; then
    echo "megahit.sh script created successfully."
else
    echo "Failed to create megahit.sh script."
    exit 1
fi

# Make the script executable
chmod +x "$MEGAHIT_ASSEMBLY_DIR/megahit.sh"
echo "Made megahit.sh executable."

# Create spades.sh script
echo "Creating spades.sh script..."
cat << EOF > "$SPADES_ASSEMBLY_DIR/spades.sh"
#!/usr/bin/bash
cd $SPADES_ASSEMBLY_DIR
READ1=\$1
READ2=\$2
OUTPUT=\$3
module load spades
mkdir -p \$OUTPUT  # Ensure the output directory exists
if [ \$? -ne 0 ]; then
    echo "Failed to create output directory \$OUTPUT"
    exit 1
fi
spades.py --meta -1 \${READ1} -2 \${READ2} -o \${OUTPUT} -t $THREADS -m $MEMORY
if [ \$? -ne 0 ]; then
    echo "SPAdes failed for \$OUTPUT"
    exit 1
fi
EOF

# Check if spades.sh was created successfully
if [[ -f "$SPADES_ASSEMBLY_DIR/spades.sh" ]]; then
    echo "spades.sh script created successfully."
else
    echo "Failed to create spades.sh script."
    exit 1
fi

# Make the script executable
chmod +x "$SPADES_ASSEMBLY_DIR/spades.sh"
echo "Made spades.sh executable."

# Create lists of forward and reverse reads
echo "Creating lists of forward and reverse reads..."
cd $BASE_DIR
ls -d $READS_DIR/*_1.fastq > "$BASE_DIR/READ1.list"
ls -d $READS_DIR/*_2.fastq > "$BASE_DIR/READ2.list"

# Check if READ1.list and READ2.list were created successfully
if [[ -f "$BASE_DIR/READ1.list" && -f "$BASE_DIR/READ2.list" ]]; then
    echo "READ1.list and READ2.list created successfully."
else
    echo "Failed to create READ1.list and/or READ2.list."
    exit 1
fi

# Read lists into arrays
echo "Reading lists into arrays..."
mapfile -t READ1 < "$BASE_DIR/READ1.list"
mapfile -t READ2 < "$BASE_DIR/READ2.list"

# Generate output directory names for Megahit
echo "Generating output directory names for Megahit..."
> "$MEGAHIT_ASSEMBLY_DIR/MEGAHIT_OUT.list"
for read1 in "${READ1[@]}"; do
    out_dir=$(basename "$read1" | sed 's/_1.fastq/_megahit_out/')
    echo "$MEGAHIT_ASSEMBLY_DIR/$out_dir" >> "$MEGAHIT_ASSEMBLY_DIR/MEGAHIT_OUT.list"
done
mapfile -t MEGAHIT_OUT < "$MEGAHIT_ASSEMBLY_DIR/MEGAHIT_OUT.list"

# Generate output directory names for Spades
echo "Generating output directory names for Spades..."
> "$SPADES_ASSEMBLY_DIR/SPADES_OUT.list"
for read1 in "${READ1[@]}"; do
    out_dir=$(basename "$read1" | sed 's/_1.fastq/_spades_out/')
    echo "$SPADES_ASSEMBLY_DIR/$out_dir" >> "$SPADES_ASSEMBLY_DIR/SPADES_OUT.list"
done
mapfile -t SPADES_OUT < "$SPADES_ASSEMBLY_DIR/SPADES_OUT.list"

# Clean up temporary files
echo "Cleaning up temporary files..."
rm "$MEGAHIT_ASSEMBLY_DIR/MEGAHIT_OUT.list" "$SPADES_ASSEMBLY_DIR/SPADES_OUT.list"

# Check if arrays are populated
if [[ ${#READ1[@]} -eq 0 || ${#READ2[@]} -eq 0 || ${#MEGAHIT_OUT[@]} -eq 0 || ${#SPADES_OUT[@]} -eq 0 ]]; then
    echo "Error: Arrays are not populated correctly."
    exit 1
else
    echo "Arrays properly populated."
fi

# Create the swarm files
echo "Creating the Megahit swarm file..."
total=${#READ1[*]}
for (( i=0; i<$total; i++ )); do
    echo "bash $MEGAHIT_ASSEMBLY_DIR/megahit.sh ${READ1[i]} ${READ2[i]} ${MEGAHIT_OUT[i]}" >> $MEGAHIT_SWARM_FILE
done

echo "Creating the Spades swarm file..."
for (( i=0; i<$total; i++ )); do
    echo "bash $SPADES_ASSEMBLY_DIR/spades.sh ${READ1[i]} ${READ2[i]} ${SPADES_OUT[i]}" >> $SPADES_SWARM_FILE
done

# Check if swarm files were created successfully
if [[ -f "$MEGAHIT_SWARM_FILE" ]]; then
    echo "Megahit swarm file created at $MEGAHIT_SWARM_FILE"
else
    echo "Failed to create Megahit swarm file."
    exit 1
fi

if [[ -f "$SPADES_SWARM_FILE" ]]; then
    echo "Spades swarm file created at $SPADES_SWARM_FILE"
else
    echo "Failed to create Spades swarm file."
    exit 1
fi

# Submit swarm jobs with user-specified settings
echo "Submitting Megahit swarm jobs with user-specified settings..."
megahit_swarmid=$(swarm -f $MEGAHIT_SWARM_FILE -g $MEMORY -t $THREADS --time $TIME --job-name "${JOB_NAME}_megahit" --logdir $MEGAHIT_LOG_DIR --module megahit)

echo "Submitting Spades swarm jobs with user-specified settings..."
spades_swarmid=$(swarm -f $SPADES_SWARM_FILE -g $MEMORY -t $THREADS --time $TIME --job-name "${JOB_NAME}_spades" --logdir $SPADES_LOG_DIR --module spades)

# Verify swarm submission
if [[ -n "$megahit_swarmid" ]]; then
    echo "Megahit swarm job submitted successfully with ID $megahit_swarmid."
else
    echo "Failed to submit Megahit swarm job."
fi

if [[ -n "$spades_swarmid" ]]; then
    echo "Spades swarm job submitted successfully with ID $spades_swarmid."
else
    echo "Failed to submit Spades swarm job."
fi
