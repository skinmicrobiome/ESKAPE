#!/bin/bash

#note: to make this work, you must use the script rename_fasta_headers.sh on your mags
#so the contigs are associated to a bin
#you must also de-replicate your eukaryotic and prokaryotic bins and concatentate them into a single fasta file for reference_genome
# Define paths and variables
forward_reads_dir="/data/proctordm/batch3/09_mapping_missing5000/test/READ1.list"
reverse_reads_dir="/data/proctordm/batch3/09_mapping_missing5000/test/READ2.list"
reference_genome="/data/proctordm/chatGPT/06_mapping/eukBacBins.derep_contam5.fa"
output_dir="/data/proctordm/chatGPT/06_mapping"

# Function to create map_and_process.sh script
#note that the script within this function was from Ryan Blaustein 
function create_map_and_process_script() {
    cat > "${output_dir}/map_and_process.sh" <<EOF
#!/bin/bash

# Arguments
reads1=\$1
reads2=\$2
ref=\$3
outprefix=\$4

# Load necessary modules
module load bwa
module load samtools

# Function to check if BWA index exists, otherwise create it
function check_bwa_index() {
    local ref=\$1
    local ref_prefix=\$(basename "\$ref" .fa)
    if [ ! -e "\${ref_prefix}.bwt" ]; then
        echo "Indexing reference genome with BWA..."
        bwa index "\$ref"
    else
        echo "BWA index already exists."
    fi
}

# Check and create BWA index if necessary
check_bwa_index "\$ref"

# Extract prefix from reference filename
ref_prefix=\$(basename "\$ref" .fa)

# Perform initial mapping
bwa mem -t 8 "\$ref" "\$reads1" "\$reads2" | samtools view -@ 7 -uS - | samtools sort -@ 7 -o "\${outprefix}_\${ref_prefix}_sorted.bam"

# Index sorted BAM file
samtools index -@ 7 "\${outprefix}_\${ref_prefix}_sorted.bam"

# Extract unique counts
samtools view -@ 7 -q 1 -f 2 -u "\${outprefix}_\${ref_prefix}_sorted.bam" -o "\${outprefix}_\${ref_prefix}_unique_sorted.bam"
samtools index -@ 7 "\${outprefix}_\${ref_prefix}_unique_sorted.bam"
samtools idxstats "\${outprefix}_\${ref_prefix}_unique_sorted.bam" > "\${outprefix}_\${ref_prefix}_unique_depth.tab"
samtools depth "\${outprefix}_\${ref_prefix}_unique_sorted.bam" > "\${outprefix}_\${ref_prefix}_unique_depth-pos.tab"
python3 parse_bwa-depth.py "\${outprefix}_\${ref_prefix}_unique_depth.tab" "\${outprefix}_\${ref_prefix}_unique_depth-pos.tab" > "\${outprefix}_\${ref_prefix}_unique.tab"
rm -rf "\${outprefix}_\${ref_prefix}_unique_sorted.bam" "\${outprefix}_\${ref_prefix}_unique_depth*"

# Extract total counts
samtools index -@ 7 "\${outprefix}_\${ref_prefix}_sorted.bam"
samtools idxstats "\${outprefix}_\${ref_prefix}_sorted.bam" > "\${outprefix}_\${ref_prefix}_depth.tab"
samtools depth "\${outprefix}_\${ref_prefix}_sorted.bam" > "\${outprefix}_\${ref_prefix}_depth-pos.tab"
python3 parse_bwa-depth.py "\${outprefix}_\${ref_prefix}_depth.tab" "\${outprefix}_\${ref_prefix}_depth-pos.tab" > "\${outprefix}_\${ref_prefix}_total.tab"
rm -rf "\${outprefix}_\${ref_prefix}_sorted.bam" "\${outprefix}_\${ref_prefix}_sorted.bam.bai" "\${outprefix}_\${ref_prefix}_depth*"

echo "Mapping and processing for \$outprefix completed."
EOF

    chmod +x "${output_dir}/map_and_process.sh"
}

# Function to create Swarm file
function generate_swarm() {
    local forward_reads="$1"
    local reverse_reads="$2"
    local reference="$3"
    local output="$4"
    
    # Read lists
    readarray -t read1_list < "$forward_reads"
    readarray -t read2_list < "$reverse_reads"
    
    # Ensure read1 and read2 lists have the same number of entries
    if [ "${#read1_list[@]}" -ne "${#read2_list[@]}" ]; then
        echo "Error: Number of entries in $forward_reads and $reverse_reads do not match."
        exit 1
    fi
    
    # Create map_and_process.sh script
    create_map_and_process_script
    
    # Create Swarm commands in output directory
    echo "Generating Swarm file..."
    rm -f "${output}/swarm_commands.swarm"  # Remove existing Swarm file if present
    
    for (( i=0; i<"${#read1_list[@]}"; i++ )); do
        read1="${read1_list[$i]}"
        read2="${read2_list[$i]}"
        base=$(basename "$read1" .fastq)
        outprefix="${output}/${base}_mapping_out"  # Adjusted output directory name
        
        # Generate Swarm command
        echo "bash ${output}/map_and_process.sh $read1 $read2 $reference $outprefix" >> "${output}/swarm_commands.swarm"
    done
    
    echo "Swarm file 'swarm_commands.swarm' generated in $output."
}

# Load necessary modules
module load python  # Assuming this is where Python 3 is available

# Generate Swarm file
generate_swarm "$forward_reads_dir" "$reverse_reads_dir" "$reference_genome" "$output_dir"

echo "Script execution completed."

