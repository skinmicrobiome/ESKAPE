#!/bin/bash

# Define directories
BINNING_DIR="/data/proctordm/ESKAPE/data/02_binning"
CONTAINER_DIR="/data/proctordm/ESKAPE/container"
OUTPUT_DIR="/data/proctordm/ESKAPE/data/03_refinement"
LOG_DIR="$OUTPUT_DIR/logs"
SIF_FILE="$CONTAINER_DIR/MAG_wf_containers_metawrap.sif"
REFINE_SCRIPT="$OUTPUT_DIR/refinement.sh"
MEGAHIT_OUTPUT_DIR="$OUTPUT_DIR/03_megahit"
SPADES_OUTPUT_DIR="$OUTPUT_DIR/03_spades"
MEGAHIT_SWARM_FILE="$MEGAHIT_OUTPUT_DIR/megahit_refinement.swarm"
SPADES_SWARM_FILE="$SPADES_OUTPUT_DIR/spades_refinement.swarm"
FAILED_SAMPLES_FILE="$OUTPUT_DIR/failed_samples.list"

# Define refinement parameters
THREADS=32
MEMORY=128
TIME="24:00:00"
MIN_COMPLETION=50
MAX_CONTAMINATION=10

# Define Swarm parameters
SWARM_JOB_NAME="refinement"
SWARM_THREADS=32
SWARM_MEMORY=128
SWARM_TIME="24:00:00"

# Ensure output and log directories exist
mkdir -p "$MEGAHIT_OUTPUT_DIR"
mkdir -p "$SPADES_OUTPUT_DIR"
mkdir -p "$LOG_DIR"

# Check if the Singularity image file exists
if [ ! -f "$SIF_FILE" ]; then
    echo "Error: Singularity image file $SIF_FILE does not exist."
    exit 1
fi

# Create refinement script
cat << EOF > "$REFINE_SCRIPT"
#!/usr/bin/bash
REFINE_OUT=\$1
metabat2_bins=\$2
maxbin2_bins=\$3
concoct_bins=\$4
SIF_FILE=\$5

module load singularity
source /usr/local/current/singularity/app_conf/sing_binds
mkdir -p \${REFINE_OUT} && cd \${REFINE_OUT}

/usr/bin/time -v singularity run \${SIF_FILE} metawrap bin_refinement -t $THREADS -m $MEMORY -o \${REFINE_OUT} -A \${metabat2_bins} -B \${maxbin2_bins} -C \${concoct_bins} -c $MIN_COMPLETION -x $MAX_CONTAMINATION
status=\$?

# Check the exit status and log the result
if [ \$status -eq 0 ]; then
    echo "Refinement job for \${REFINE_OUT} finished successfully." >> \${REFINE_OUT}/refinement.log
else
    echo "Refinement job for \${REFINE_OUT} failed with status \$status." >> \${REFINE_OUT}/refinement.log
fi

# Log resource usage
echo "Resource usage for refinement job \${REFINE_OUT}:" >> \${REFINE_OUT}/resource_usage.log
grep "Elapsed (wall clock) time" \${REFINE_OUT}/time.log >> \${REFINE_OUT}/resource_usage.log
grep "Maximum resident set size" \${REFINE_OUT}/time.log >> \${REFINE_OUT}/resource_usage.log
grep "User time (seconds)" \${REFINE_OUT}/time.log >> \${REFINE_OUT}/resource_usage.log
grep "System time (seconds)" \${REFINE_OUT}/time.log >> \${REFINE_OUT}/resource_usage.log
EOF

# Make the script executable
chmod +x "$REFINE_SCRIPT"

# Clear previous lists
> "$MEGAHIT_OUTPUT_DIR/megahit_metabat2_bins.list"
> "$MEGAHIT_OUTPUT_DIR/megahit_maxbin2_bins.list"
> "$MEGAHIT_OUTPUT_DIR/megahit_concoct_bins.list"
> "$SPADES_OUTPUT_DIR/spades_metabat2_bins.list"
> "$SPADES_OUTPUT_DIR/spades_maxbin2_bins.list"
> "$SPADES_OUTPUT_DIR/spades_concoct_bins.list"
> "$MEGAHIT_OUTPUT_DIR/megahit_metabat2_samples.list"
> "$MEGAHIT_OUTPUT_DIR/megahit_maxbin2_samples.list"
> "$MEGAHIT_OUTPUT_DIR/megahit_concoct_samples.list"
> "$SPADES_OUTPUT_DIR/spades_metabat2_samples.list"
> "$SPADES_OUTPUT_DIR/spades_maxbin2_samples.list"
> "$SPADES_OUTPUT_DIR/spades_concoct_samples.list"
> "$MEGAHIT_OUTPUT_DIR/megahit_common_samples.list"
> "$SPADES_OUTPUT_DIR/spades_common_samples.list"
> "$MEGAHIT_OUTPUT_DIR/megahit_failed_samples.list"
> "$SPADES_OUTPUT_DIR/spades_failed_samples.list"
> "$MEGAHIT_OUTPUT_DIR/megahit_REFINEMENT_OUT.list"
> "$SPADES_OUTPUT_DIR/spades_REFINEMENT_OUT.list"

# List bin directories for MegaHit and SPAdes
ls -d $BINNING_DIR/megahit/*/metabat2_bins > "$MEGAHIT_OUTPUT_DIR/megahit_metabat2_bins.list"
ls -d $BINNING_DIR/megahit/*/maxbin2_bins > "$MEGAHIT_OUTPUT_DIR/megahit_maxbin2_bins.list"
ls -d $BINNING_DIR/megahit/*/concoct_bins > "$MEGAHIT_OUTPUT_DIR/megahit_concoct_bins.list"

ls -d $BINNING_DIR/spades/*/metabat2_bins > "$SPADES_OUTPUT_DIR/spades_metabat2_bins.list"
ls -d $BINNING_DIR/spades/*/maxbin2_bins > "$SPADES_OUTPUT_DIR/spades_maxbin2_bins.list"
ls -d $BINNING_DIR/spades/*/concoct_bins > "$SPADES_OUTPUT_DIR/spades_concoct_bins.list"

# Extract sample names from paths
while read -r path; do
    basename "$(dirname "$path")" | sed 's/_binning_out//' >> "$MEGAHIT_OUTPUT_DIR/megahit_metabat2_samples.list"
done < "$MEGAHIT_OUTPUT_DIR/megahit_metabat2_bins.list"

while read -r path; do
    basename "$(dirname "$path")" | sed 's/_binning_out//' >> "$MEGAHIT_OUTPUT_DIR/megahit_maxbin2_samples.list"
done < "$MEGAHIT_OUTPUT_DIR/megahit_maxbin2_bins.list"

while read -r path; do
    basename "$(dirname "$path")" | sed 's/_binning_out//' >> "$MEGAHIT_OUTPUT_DIR/megahit_concoct_samples.list"
done < "$MEGAHIT_OUTPUT_DIR/megahit_concoct_bins.list"

while read -r path; do
    basename "$(dirname "$path")" | sed 's/_binning_out//' >> "$SPADES_OUTPUT_DIR/spades_metabat2_samples.list"
done < "$SPADES_OUTPUT_DIR/spades_metabat2_bins.list"

while read -r path; do
    basename "$(dirname "$path")" | sed 's/_binning_out//' >> "$SPADES_OUTPUT_DIR/spades_maxbin2_samples.list"
done < "$SPADES_OUTPUT_DIR/spades_maxbin2_bins.list"

while read -r path; do
    basename "$(dirname "$path")" | sed 's/_binning_out//' >> "$SPADES_OUTPUT_DIR/spades_concoct_samples.list"
done < "$SPADES_OUTPUT_DIR/spades_concoct_bins.list"

# Find common samples for MegaHit
comm -12 <(sort "$MEGAHIT_OUTPUT_DIR/megahit_metabat2_samples.list") <(sort "$MEGAHIT_OUTPUT_DIR/megahit_maxbin2_samples.list") > "$MEGAHIT_OUTPUT_DIR/megahit_common_samples.temp"
comm -12 <(sort "$MEGAHIT_OUTPUT_DIR/megahit_common_samples.temp") <(sort "$MEGAHIT_OUTPUT_DIR/megahit_concoct_samples.list") > "$MEGAHIT_OUTPUT_DIR/megahit_common_samples.list"

# Find common samples for SPAdes
comm -12 <(sort "$SPADES_OUTPUT_DIR/spades_metabat2_samples.list") <(sort "$SPADES_OUTPUT_DIR/spades_maxbin2_samples.list") > "$SPADES_OUTPUT_DIR/spades_common_samples.temp"
comm -12 <(sort "$SPADES_OUTPUT_DIR/spades_common_samples.temp") <(sort "$SPADES_OUTPUT_DIR/spades_concoct_samples.list") > "$SPADES_OUTPUT_DIR/spades_common_samples.list"

# Remove temp files
rm "$MEGAHIT_OUTPUT_DIR/megahit_common_samples.temp"
rm "$SPADES_OUTPUT_DIR/spades_common_samples.temp"

# Identify failed samples and report them
comm -23 <(sort "$MEGAHIT_OUTPUT_DIR/megahit_metabat2_samples.list") <(sort "$MEGAHIT_OUTPUT_DIR/megahit_common_samples.list") > "$MEGAHIT_OUTPUT_DIR/megahit_failed_samples.list"
comm -23 <(sort "$SPADES_OUTPUT_DIR/spades_metabat2_samples.list") <(sort "$SPADES_OUTPUT_DIR/spades_common_samples.list") > "$SPADES_OUTPUT_DIR/spades_failed_samples.list"

# Combine failed samples
cat "$MEGAHIT_OUTPUT_DIR/megahit_failed_samples.list" "$SPADES_OUTPUT_DIR/spades_failed_samples.list" | sort | uniq > "$FAILED_SAMPLES_FILE"

# Check if there are failed samples and report
if [[ -s "$FAILED_SAMPLES_FILE" ]]; then
    echo "Failed samples detected:"
    cat "$FAILED_SAMPLES_FILE"
else
    echo "No failed samples."
fi

# Generate output directory names
while read -r sample; do
    echo "$MEGAHIT_OUTPUT_DIR/${sample}_megahit_refinement_out" >> "$MEGAHIT_OUTPUT_DIR/megahit_REFINEMENT_OUT.list"
done < "$MEGAHIT_OUTPUT_DIR/megahit_common_samples.list"

while read -r sample; do
    echo "$SPADES_OUTPUT_DIR/${sample}_spades_refinement_out" >> "$SPADES_OUTPUT_DIR/spades_REFINEMENT_OUT.list"
done < "$SPADES_OUTPUT_DIR/spades_common_samples.list"

# Generate Swarm files for MegaHit and SPAdes from common samples list
mapfile -t megahit_common_samples < "$MEGAHIT_OUTPUT_DIR/megahit_common_samples.list"
mapfile -t spades_common_samples < "$SPADES_OUTPUT_DIR/spades_common_samples.list"

# Create Swarm file for MegaHit
> "$MEGAHIT_SWARM_FILE"
total=${#megahit_common_samples[*]}
for (( i=0; i<$total; i++ )); do
    sample="${megahit_common_samples[i]}"
    output_dir="$MEGAHIT_OUTPUT_DIR/${sample}_megahit_refinement_out"
    metabat2_bins="$BINNING_DIR/megahit/${sample}_binning_out/metabat2_bins"
    maxbin2_bins="$BINNING_DIR/megahit/${sample}_binning_out/maxbin2_bins"
    concoct_bins="$BINNING_DIR/megahit/${sample}_binning_out/concoct_bins"
    echo "bash $REFINE_SCRIPT $output_dir $metabat2_bins $maxbin2_bins $concoct_bins $SIF_FILE" >> "$MEGAHIT_SWARM_FILE"
done

# Create Swarm file for SPAdes
> "$SPADES_SWARM_FILE"
total=${#spades_common_samples[*]}
for (( i=0; i<$total; i++ )); do
    sample="${spades_common_samples[i]}"
    output_dir="$SPADES_OUTPUT_DIR/${sample}_spades_refinement_out"
    metabat2_bins="$BINNING_DIR/spades/${sample}_binning_out/metabat2_bins"
    maxbin2_bins="$BINNING_DIR/spades/${sample}_binning_out/maxbin2_bins"
    concoct_bins="$BINNING_DIR/spades/${sample}_binning_out/concoct_bins"
    echo "bash $REFINE_SCRIPT $output_dir $metabat2_bins $maxbin2_bins $concoct_bins $SIF_FILE" >> "$SPADES_SWARM_FILE"
done

# Check if Swarm files were created successfully
if [[ -f "$MEGAHIT_SWARM_FILE" ]]; then
    echo "Megahit Swarm file created at $MEGAHIT_SWARM_FILE"
else
    echo "Failed to create Megahit Swarm file."
    exit 1
fi

if [[ -f "$SPADES_SWARM_FILE" ]]; then
    echo "Spades Swarm file created at $SPADES_SWARM_FILE"
else
    echo "Failed to create Spades Swarm file."
    exit 1
fi

# Submit the Megahit Swarm job
echo "Submitting the Megahit Swarm job..."
megahit_swarmid=$(swarm -f $MEGAHIT_SWARM_FILE --job-name ${SWARM_JOB_NAME}_megahit -t $SWARM_THREADS -g $SWARM_MEMORY --time $SWARM_TIME --logdir $LOG_DIR)

# Submit the Spades Swarm job
echo "Submitting the Spades Swarm job..."
spades_swarmid=$(swarm -f $SPADES_SWARM_FILE --job-name ${SWARM_JOB_NAME}_spades -t $SWARM_THREADS -g $SWARM_MEMORY --time $SWARM_TIME --logdir $LOG_DIR)

# Check if Swarm jobs were submitted successfully
if [[ -z "$megahit_swarmid" ]]; then
    echo "Failed to submit Megahit Swarm jobs."
    exit 1
else
    echo "Megahit Swarm jobs submitted with ID $megahit_swarmid"
fi

if [[ -z "$spades_swarmid" ]]; then
    echo "Failed to submit Spades Swarm jobs."
    exit 1
else
    echo "Spades Swarm jobs submitted with ID $spades_swarmid"
fi

# Wait for Megahit Swarm jobs to finish
echo "Waiting for Megahit Swarm jobs to finish..."
while squeue -u $USER | grep -q $megahit_swarmid; do sleep 10; done

# Wait for Spades Swarm jobs to finish
echo "Waiting for Spades Swarm jobs to finish..."
while squeue -u $USER | grep -q $spades_swarmid; do sleep 10; done

echo "Refinement process completed."

# Analyze resource usage logs
echo "Analyzing resource usage logs..."
grep "Elapsed (wall clock) time" $MEGAHIT_OUTPUT_DIR/*/resource_usage.log
grep "Maximum resident set size" $MEGAHIT_OUTPUT_DIR/*/resource_usage.log
grep "User time (seconds)" $MEGAHIT_OUTPUT_DIR/*/resource_usage.log
grep "System time (seconds)" $MEGAHIT_OUTPUT_DIR/*/resource_usage.log

grep "Elapsed (wall clock) time" $SPADES_OUTPUT_DIR/*/resource_usage.log
grep "Maximum resident set size" $SPADES_OUTPUT_DIR/*/resource_usage.log
grep "User time (seconds)" $SPADES_OUTPUT_DIR/*/resource_usage.log
grep "System time (seconds)" $SPADES_OUTPUT_DIR/*/resource_usage.log
