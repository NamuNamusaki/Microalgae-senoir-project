#!/bin/bash

# --- 1. Configuration ---
# Update this to the specific species folder (e.g., Haematococcus_pluvialis)
SPECIES_NAME="Chlorella_sp" 

# Base paths for your project on sbi_cold
DATA_DIR="/home_sbi_cold/thana.yens/senior_project/transcriptome_data/$SPECIES_NAME"
INDEX_DIR="/home_sbi_cold/thana.yens/senior_project/result_all/genomic/rsem_index/$SPECIES_NAME"
STAR_PATH="/home_sbi_cold/thana.yens/STAR-2.7.11b/source"
THREADS=16

# The specific RSEM index prefix
RSEM_INDEX="$INDEX_DIR/${SPECIES_NAME}_rep"

# --- 2. The Calculation Loop ---
for sample_dir in "$DATA_DIR"/SRR*/; do
    
    # Extract the Accession ID (e.g., SRR29906374)
    ACCESSION=$(basename "$sample_dir")
    
    echo "======================================================="
    echo "Processing Accession: $ACCESSION"
    echo "======================================================="

    # Define the potential FASTQ files
    R1="$sample_dir/${ACCESSION}_1.fastq"
    R2="$sample_dir/${ACCESSION}_2.fastq"
    SINGLE="$sample_dir/${ACCESSION}.fastq"
    OUT_PREFIX="$sample_dir/${ACCESSION}_rsem"

    # 1. Check if calculation is already done
    if [ -f "${OUT_PREFIX}.genes.results" ]; then
        echo "   [Skip] Results already exist for $ACCESSION."
        continue
    fi

    # 2. Cleanup any failed temp folders from previous crashes
    rm -rf "${OUT_PREFIX}.temp"

    # --- 3. Detection & Execution Logic ---

    # CASE A: Paired-End Detection
    if [[ -f "$R1" && -f "$R2" ]]; then
        echo "   [Action] Running Paired-End RSEM quantification..."
        
        rsem-calculate-expression \
            --star \
            --star-path "$STAR_PATH" \
            --paired-end \
            --num-threads "$THREADS" \
            --no-bam-output \
            --estimate-rspd \
            --append-names \
            "$R1" "$R2" \
            "$RSEM_INDEX" \
            "$OUT_PREFIX"

    # CASE B: Single-End Detection (Accessions ending in .fastq)
    elif [[ -f "$SINGLE" ]]; then
        echo "   [Action] Running Single-End RSEM quantification..."
        
        rsem-calculate-expression \
            --star \
            --star-path "$STAR_PATH" \
            --num-threads "$THREADS" \
            --no-bam-output \
            --estimate-rspd \
            --append-names \
            "$SINGLE" \
            "$RSEM_INDEX" \
            "$OUT_PREFIX"

    # CASE C: Files Missing
    else
        echo "   [Error] No valid FASTQ files found for $ACCESSION in $sample_dir"
        continue
    fi

    # Check for success
    if [ $? -eq 0 ]; then
        echo "   [Success] TPM calculation complete for $ACCESSION."
    else
        echo "   [Error] RSEM failed for $ACCESSION. Check log files."
    fi
done

echo "All calculations for $SPECIES_NAME are finished!"
