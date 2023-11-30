#!/bin/bash

# Load necessary modules
module load fastqc
module load cutadapt

# Directory for FastQC results
mkdir -p ./fastqc_results/

# FastQC for quality control
# Generates a report to select appropriate adapters and primers for trimming
fastqc -o ./fastqc_results/ *.fastq

# Adapter trimming with Cutadapt
# Adapters and primers were selected based on FastQC report analysis
for file in *.fastq; do
    cutadapt -a AGATCGGAAGAG -a CGCCTTGGCCGT -g GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT \
    --minimum-length 1 -o "trimmed_${file}" "${file}"
done

echo "Quality control and trimming complete based on FastQC report analysis."
