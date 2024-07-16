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
LOG_DIR="$BASE_DIR/data/01_assembly/swarm_logs"

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
mkdir -p "$MEGAHIT_ASSEMBLY_DIR"
mkdir -p "$SPADES_ASSEMBLY_DIR"
mkdir -p "$LOG_DIR"

# Create megahit.sh script
echo "Creating megahit.sh script..."
cat << EOF > "$MEGAHIT_ASSEMBLY_DIR/megahit.sh"
#!/usr/bin/bash
cd $MEGAHIT_ASSEMBLY_DIR
READ1=\$1
READ2=\$2
OUTPUT=\$3
module load megahit
/usr/bin/time -v megahit -1 \${READ1} -2 \${READ2} -o \${OUTPUT} --memory $MEMORY -t $THREADS 2> \${OUTPUT}/time.log
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
/usr/bin/time -v spades.py --meta -1 \${READ1} -2 \${READ2} -o \${OUTPUT} -t $THREADS -m $MEMORY 2> \${OUTPUT}/time.log
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
megahit_swarmid=$(swarm -f $MEGAHIT_SWARM_FILE -g $MEMORY -t $THREADS --time $TIME --job-name "${JOB_NAME}_megahit" --logdir $LOG_DIR)

echo "Submitting Spades swarm jobs with user-specified settings..."
spades_swarmid=$(swarm -f $SPADES_SWARM_FILE -g $MEMORY -t $THREADS --time $TIME --job-name "${JOB_NAME}_spades" --logdir $LOG_DIR)

# Check if swarm jobs were submitted successfully
if [[ -z "$megahit_swarmid" ]]; then
    echo "Failed to submit Megahit swarm jobs."
    exit 1
else
    echo "Megahit swarm jobs submitted with ID $megahit_swarmid"
fi

if [[ -z "$spades_swarmid" ]]; then
    echo "Failed to submit Spades swarm jobs."
    exit 1
else
    echo "Spades swarm jobs submitted with ID $spades_swarmid"
fi

# Wait for swarm jobs to finish
echo "Waiting for Megahit swarm jobs to finish..."
while squeue -u $USER -j $megahit_swarmid > /dev/null 2>&1; do
    sleep 60
done 

echo "Waiting for Spades swarm jobs to finish..."
while squeue -u $USER -j $spades_swarmid > /dev/null 2>&1; do
    sleep 60
done

# Check if the assembly directories exist before finding final contigs
if [[ -d "$MEGAHIT_ASSEMBLY_DIR" ]]; then
    # Create list of final contigs for Megahit
    echo "Creating list of final contigs for Megahit..."
    find $MEGAHIT_ASSEMBLY_DIR -name "final.contigs.fa" > $MEGAHIT_FINAL_CONTIGS_LIST

    # Check if final contigs list was created successfully
    if [[ -s "$MEGAHIT_FINAL_CONTIGS_LIST" ]]; then
        echo "Final contigs list created at $MEGAHIT_FINAL_CONTIGS_LIST"
    else
        echo "No final contigs files found for Megahit."
    fi
else
    echo "Megahit assembly directory $MEGAHIT_ASSEMBLY_DIR does not exist."
    exit 1
fi

if [[ -d "$SPADES_ASSEMBLY_DIR" ]]; then
    # Create list of final contigs for Spades
    echo "Creating list of final contigs for Spades..."
    find $SPADES_ASSEMBLY_DIR -name "contigs.fasta" > $SPADES_FINAL_CONTIGS_LIST

    # Check if final contigs list was created successfully
    if [[ -s "$SPADES_FINAL_CONTIGS_LIST" ]]; then
        echo "Final contigs list created at $SPADES_FINAL_CONTIGS_LIST"
    else
        echo "No final contigs files found for Spades."
    fi
else
    echo "Spades assembly directory $SPADES_ASSEMBLY_DIR does not exist."
    exit 1
fi

# Analyze time.log files to determine resource usage
echo "Analyzing time.log files for resource usage..."
grep "Elapsed (wall clock) time" $MEGAHIT_ASSEMBLY_DIR/*/time.log
grep "Maximum resident set size" $MEGAHIT_ASSEMBLY_DIR/*/time.log
grep "Elapsed (wall clock) time" $SPADES_ASSEMBLY_DIR/*/time.log
grep "Maximum resident set size" $SPADES_ASSEMBLY_DIR/*/time.log
