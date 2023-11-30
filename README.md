# Transcriptomics_Analysis-or-COVID19_Response_Transcriptome
Project Description
This project focuses on understanding the transcriptomic changes in human respiratory cells in response to SARS-CoV-2 infection, the virus responsible for the COVID-19 pandemic. By comparing miRNA expression levels in cultivated human respiratory cells at different time points, we aim to shed light on the virus's impact on cellular functions and the immune response. Our analysis is based on data from NCBI BioProject PRJNA901149.

Assignment Objectives:
Align and quantify transcriptome profiles for control (mock) and SARS-CoV-2 infected respiratory cells at 24 and 72 hours post-infection.
Identify differentially expressed genes (DEGs) across conditions and time points.
Conduct GO Term enrichment analysis for identified DEGs.
Methodology:
The analysis workflow includes quality control of raw sequence data, trimming of adapter sequences, alignment to the human reference genome, quantification of gene expression levels, differential expression analysis, and GO Term enrichment analysis.

Tools and Resources Used:
FastQC for quality assessment of raw sequence data.
Cutadapt and Trimmomatic for trimming adapter sequences and low-quality bases.
Bowtie2 for sequence alignment to the human reference genome.
DESeq2 for differential gene expression analysis.
ClusterProfiler for GO Term enrichment analysis.
The human genome build utilized for alignment is GRCh38.

Results:
The output includes a comprehensive set of files detailing gene expression levels across all conditions and samples, alongside the results of differential expression analysis and GO Term enrichment.

Please refer to the individual script files provided in this repository for step-by-step instructions on executing the analysis pipeline.
