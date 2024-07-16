#!/bin/bash

# Set up directories and paths
BASE_DIR="/data/proctordm/chatGPT"
QC_PROKARYOTE_DIR="$BASE_DIR/data/04_QC/QC_prokaryote_all"
QC_EUKARYOTE_DIR="$BASE_DIR/data/04_QC/eukaryote_all"
QC_EUKARYOTE_MEGAHIT_DIR="$QC_EUKARYOTE_DIR/megahit"
QC_EUKARYOTE_SPADES_DIR="$QC_EUKARYOTE_DIR/spades"
MASH_OUTPUT_DIR="$BASE_DIR/data/05_euk"
FUNGI_GENBANK_DIR="/data/$USER/metagenomes_mags/data/databases/fungi_genbank"
FUNGI_MASH_SKETCH="/data/$USER/metagenomes_mags/data/databases/fungi.genbank.msh"
BESTMASHHITS_OUTPUT_DIR="$BASE_DIR/data/05_euk"
MUMMER_DIR="$BASE_DIR/data/05_euk/mummer"
GTDBTK_MEGAHIT_OUTPUT="$QC_PROKARYOTE_DIR/qc_megahit/gtdbk"
GTDBTK_SPADES_OUTPUT="$QC_PROKARYOTE_DIR/qc_spades/gtdbk"
PROKARYOTE_MEGAHIT_INPUT="$QC_PROKARYOTE_DIR/qc_megahit/input"
PROKARYOTE_SPADES_INPUT="$QC_PROKARYOTE_DIR/qc_spades/input"
cpus_per_task="${SLURM_CPUS_PER_TASK:-1}"  # Default to 1 if not set

# Ensure necessary directories exist
mkdir -p "$BASE_DIR/data" "$QC_PROKARYOTE_DIR" "$QC_EUKARYOTE_DIR" "$QC_EUKARYOTE_MEGAHIT_DIR" "$QC_EUKARYOTE_SPADES_DIR" "$MASH_OUTPUT_DIR" "$FUNGI_GENBANK_DIR" "$BESTMASHHITS_OUTPUT_DIR" "$MUMMER_DIR" "$GTDBTK_MEGAHIT_OUTPUT" "$GTDBTK_SPADES_OUTPUT"

# Determine memory per task using Slurm environment
if [ -n "$SLURM_MEM_PER_CPU" ]; then
    mem=$(echo "$SLURM_MEM_PER_CPU * $cpus_per_task" | bc)m
elif [ -n "$SLURM_MEM_PER_NODE" ]; then
    mem="$SLURM_MEM_PER_NODE"
else
    mem="128g"  # Default value if Slurm environment variables are not available
fi

time="36:00:00"

# Function to log errors
log_error() {
    local msg="$1"
    echo "ERROR: $msg" >&2
}

# Function to set up fungal GenBank database and create MASH sketch if not already done
setup_fungi_genbank() {
    mkdir -p "$FUNGI_GENBANK_DIR"
    cd "$FUNGI_GENBANK_DIR"
    
    # Download the fungal GenBank database
    ncbi-genome-download fungi --formats fasta --section genbank
    
    cd "$FUNGI_GENBANK_DIR/genbank/fungi"
    mkdir -p fastas
    cp **/*fna.gz fastas/
    gunzip fastas/*.gz
    
    # Include species specific to your dataset
    cp /data/Segrelab/data/all_genomes/Caur_296.Fg_acrtq.spades.pilon.fasta fastas/
    cp /data/Segrelab/data/all_genomes/Caur_2117.An_acuzh.spades.pilon.fasta fastas/
    cp /data/Segrelab/data/all_genomes/Caur_283.Ic_acrty.spades.pilon.fasta fastas/
    
    # Include SMGC fungal species
    cp /data/Segrelab/data/zoo/SMGCe/*fa fastas/
    
    # Create the MASH sketch
    module load mash
    cd /data/$USER/metagenomes_mags/data/databases
    mash sketch -o "$FUNGI_MASH_SKETCH" "$FUNGI_GENBANK_DIR/genbank/fungi/fastas/*"
}

# Function to run GTDB-Tk on prokaryotes
run_gtdbtk_prokaryotes() {
    source myconda
    conda activate gtdbtk-2.1.1

    # Run GTDB-Tk classification on MegaHit assemblies
    mkdir -p "$GTDBTK_MEGAHIT_OUTPUT"
    echo "Running GTDB-Tk on MegaHit assemblies..."
    echo "Input directory: $PROKARYOTE_MEGAHIT_INPUT"
    echo "Output directory: $GTDBTK_MEGAHIT_OUTPUT"
    gtdbtk classify_wf --cpus "$cpus_per_task" --genome_dir "$PROKARYOTE_MEGAHIT_INPUT" --out_dir "$GTDBTK_MEGAHIT_OUTPUT" -x fa
    if [ $? -eq 0 ]; then
        echo "GTDB-Tk classification for MegaHit assemblies completed successfully."
    else
        log_error "GTDB-Tk classification for MegaHit assemblies failed."
    fi

    # Run GTDB-Tk classification on SPAdes assemblies
    mkdir -p "$GTDBTK_SPADES_OUTPUT"
    echo "Running GTDB-Tk on SPAdes assemblies..."
    echo "Input directory: $PROKARYOTE_SPADES_INPUT"
    echo "Output directory: $GTDBTK_SPADES_OUTPUT"
    gtdbtk classify_wf --cpus "$cpus_per_task" --genome_dir "$PROKARYOTE_SPADES_INPUT" --out_dir "$GTDBTK_SPADES_OUTPUT" -x fa
    if [ $? -eq 0 ]; then
        echo "GTDB-Tk classification for SPAdes assemblies completed successfully."
    else
        log_error "GTDB-Tk classification for SPAdes assemblies failed."
    fi
}

# Function to perform MASH sketch and comparison for eukaryotes
run_mash_eukaryotes() {
    # Ensure the MASH sketch file exists, if not create it
    if [ ! -f "$FUNGI_MASH_SKETCH" ]; then
        setup_fungi_genbank
    fi

    for assembler in megahit spades; do
        if [ "$assembler" == "megahit" ]; then
            eukaryote_input_dir="$QC_EUKARYOTE_MEGAHIT_DIR"
        elif [ "$assembler" == "spades" ]; then
            eukaryote_input_dir="$QC_EUKARYOTE_SPADES_DIR"
        fi
        mash_output="$MASH_OUTPUT_DIR/${assembler}_fungi.genbank_mash_out.txt"
        bestmashhits_output="$BESTMASHHITS_OUTPUT_DIR/${assembler}_bestmashhits.fungi.genbank.msh_manual.txt"

        # MASH sketch comparison against fungal genbank for eukaryotic bins
        module load mash
        if [ ! -f "$mash_output" ]; then
            echo "Running mash dist -p $cpus_per_task $FUNGI_MASH_SKETCH $eukaryote_input_dir/*.fa > $mash_output"
            
            # Check if input files exist
            if ls $eukaryote_input_dir/*.fa 1> /dev/null 2>&1; then
                mash dist -p "$cpus_per_task" "$FUNGI_MASH_SKETCH" "$eukaryote_input_dir"/*.fa > "$mash_output" 2> mash_error.log
                if [ $? -eq 0 ]; then
                    echo "MASH sketch comparison against fungal genbank completed for $assembler."
                else
                    log_error "MASH sketch comparison failed for $assembler. Check mash_error.log for details."
                    continue
                fi
            else
                log_error "No input files found in $eukaryote_input_dir for $assembler."
                continue
            fi
        else
            echo "MASH output file '$mash_output' already exists for $assembler. Skipping MASH comparison."
        fi

        # Check if MASH output file was created successfully and is not empty
        if [ -s "$mash_output" ]; then
            # Process MASH output to find best hits for each MAG
            awk 'BEGIN { FS="\t"; OFS="\t"; prev_MAG="" }
                 { if ($2 != prev_MAG) { print prev_line; prev_line = $0; prev_MAG = $2 }
                   else if ($4 < min_pvalue) { min_pvalue = $4; prev_line = $0 }
                 }
                 END { print prev_line }' "$mash_output" | awk '{print $1, $2, $3, $4, $5, $6, $1}' OFS="\t" > "$bestmashhits_output"

            echo "Best MASH hits file generated for $assembler: $bestmashhits_output"
        else
            log_error "MASH output file '$mash_output' not found or empty for $assembler."
        fi
    done
}

# Function to run MUMmer on eukaryote MAGs
run_mummer_eukaryotes() {
    cd "$MUMMER_DIR"
    module load mummer

    # Create MUMmer script
    cat << EOF > mummer.sh
#!/usr/bin/bash
REFERENCE=\$1
MAG=\$2
OUT=\$(basename \${MAG} .fa)
module load mummer
mkdir -p \${OUT}_mummer_output && cd \${OUT}_mummer_output
dnadiff \${REFERENCE} \${MAG} -p \${OUT}
file=\${OUT}.rev
cat *report | grep -i "^AlignedBases" > \$file
cat *report | grep -i "^AvgIdentity" | awk 'NR%2!=0' >> \$file
cp \$file $MUMMER_DIR/\${OUT}.rev
EOF

    # Create swarm file for each assembler
    for assembler in megahit spades; do
        input="$BESTMASHHITS_OUTPUT_DIR/${assembler}_bestmashhits.fungi.genbank.msh_manual.txt"
        if [ ! -f "$input" ]; then
            log_error "Input file $input not found for $assembler."
            continue
        fi

        REFERENCE=($(awk '{print $1}' $input))
        MAG=($(awk '{print $2}' $input))

        # Debugging statements to check REFERENCE and MAG values
        echo "DEBUG: Assembler: $assembler"
        echo "DEBUG: REFERENCE: ${REFERENCE[@]}"
        echo "DEBUG: MAG: ${MAG[@]}"

        total=${#REFERENCE[*]}
        for ((i=0; i<$(($total)); i++)); do
            OUT=$(basename ${MAG[i]} .fa)
            echo "bash mummer.sh ${REFERENCE[i]} ${MAG[i]} $OUT"
        done > mummer_${assembler}.swarm

        # Run swarm
        swarm -f mummer_${assembler}.swarm --job-name mummer_${assembler} -t 16 -g 98 --time 2:00:00 -b 2 --logdir mummer_swarm_out
    done
}

# Function to aggregate MUMmer reports
aggregate_mummer_reports() {
    cd "$MUMMER_DIR"

    # Use globbing to find all .rev files and store them in an array
    mummerOut=($MUMMER_DIR/**/*.rev)

    # Check if any .rev files were found
    if [ ${#mummerOut[@]} -eq 0 ]; then
        echo "No .rev files found in $MUMMER_DIR"
        return 1
    fi

    # Initialize arrays to store the data
    AlignedBases=()
    AvgIdentity=()
    mummerBin=()

    # Read each MUMmer output file and extract data
    for file in "${mummerOut[@]}"; do
        alignedBases=$(grep -i "^AlignedBases" "$file")
        avgIdentity=$(grep -i "^AvgIdentity" "$file" | awk 'NR%2!=0')
        bin=$(basename "$file" .rev)
        AlignedBases+=("$alignedBases")
        AvgIdentity+=("$avgIdentity")
        mummerBin+=("$bin")
    done

    # Create a CSV file with the aggregated data
    output_file="mummer_merged.csv"
    {
        echo "AlignedBases,AvgIdentity,Bin"
        for i in "${!AlignedBases[@]}"; do
            echo "${AlignedBases[i]},${AvgIdentity[i]},${mummerBin[i]}"
        done
    } > "$output_file"

    echo "MUMmer reports aggregation completed and saved to $output_file."
}

# Main script execution starts here
run_gtdbtk_prokaryotes
run_mash_eukaryotes
run_mummer_eukaryotes
aggregate_mummer_reports

echo "Script execution completed successfully."
#you have to merge the mash, mummer and ncbi taxonomy report to get the fungal tax file
