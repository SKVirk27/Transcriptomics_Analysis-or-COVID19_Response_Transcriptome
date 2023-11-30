#!/bin/bash

# Load Bowtie2 and Samtools modules if on a shared compute environment
module load bowtie2
module load samtools

# Bowtie2 index path - replace with the actual path after downloading the index
# The index used is the GRCh38 no-alt analysis set from NCBI, downloaded via HTTPS URLs
BT2_INDEX="/path/to/bowtie2/index/GRCh38_noalt_as"

# Base directory containing trimmed FASTQ files
BASE_DIR="/path/to/trimmed_fastq"

# Output directory for SAM and BAM files
OUTPUT_DIR="./aligned_data"
mkdir -p $OUTPUT_DIR

# Align reads and convert to sorted BAM
for file in ${BASE_DIR}/*_trimmed.fastq; do
    fname=$(basename "$file" "_trimmed.fastq")
    echo "Aligning $fname with Bowtie2..."
    bowtie2 -x $BT2_INDEX -U $file -S ${OUTPUT_DIR}/${fname}.sam
    echo "Converting SAM to sorted BAM for $fname..."
    samtools view -bS ${OUTPUT_DIR}/${fname}.sam | samtools sort -o ${OUTPUT_DIR}/${fname}_sorted.bam
    samtools index ${OUTPUT_DIR}/${fname}_sorted.bam
    rm ${OUTPUT_DIR}/${fname}.sam
done

echo "Alignment with Bowtie2 and conversion to BAM complete for all samples."
