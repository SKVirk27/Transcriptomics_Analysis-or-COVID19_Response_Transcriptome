# Transcriptomics Analysis of COVID-19 Response Transcriptome

## Project Description
This project aims to understand the transcriptomic changes in human respiratory cells in response to SARS-CoV-2 infection, the virus responsible for the COVID-19 pandemic. By comparing miRNA expression levels in cultivated human respiratory cells at different time points, we aim to shed light on the virus's impact on cellular functions and the immune response. The analysis is based on data from NCBI BioProject PRJNA901149.

## Assignment Objectives
- Align and quantify transcriptome profiles for control (mock) and SARS-CoV-2 infected respiratory cells at 24 and 72 hours post-infection.
- Identify differentially expressed genes (DEGs) across conditions and time points.
- Conduct GO Term enrichment analysis for identified DEGs.

## Methodology
The analysis workflow includes:
1. **Quality Control**: Assessing raw sequence data quality using FastQC.
2. **Trimming**: Removing adapter sequences and low-quality bases using Cutadapt and Trimmomatic.
3. **Alignment**: Aligning sequences to the human reference genome (GRCh38) using Bowtie2.
4. **Quantification**: Quantifying gene expression levels.
5. **Differential Expression Analysis**: Identifying DEGs using DESeq2.
6. **GO Term Enrichment Analysis**: Analyzing GO Terms for identified DEGs using ClusterProfiler.

## Tools and Resources Used
- **FastQC**: Quality assessment of raw sequence data.
- **Cutadapt and Trimmomatic**: Trimming adapter sequences and low-quality bases.
- **Bowtie2**: Sequence alignment to the human reference genome.
- **DESeq2**: Differential gene expression analysis.
- **ClusterProfiler**: GO Term enrichment analysis.
- **Human Genome Build GRCh38**: Used for alignment.

## Results
The output includes:
- A comprehensive set of files detailing gene expression levels across all conditions and samples.
- Results of differential expression analysis.
- GO Term enrichment results.

## Repository Structure
```
.
├── Differential_expression_using_DESeq2.R
├── GO_condition.csv
├── GO_time.csv
├── Gene_ontology.R
├── Mapping_with_Bowtie.sh
├── Quality_mapping_scores.xlsx
├── README.md
├── Readme.docx
├── Trimming.sh
├── barplot_condition.pdf
├── barplot_time.pdf
├── ego_condition.pdf
├── ego_time.pdf
```

### Key Files and Scripts
- **Differential_expression_using_DESeq2.R**: Script for differential gene expression analysis using DESeq2.
- **Gene_ontology.R**: Script for conducting GO Term enrichment analysis.
- **Mapping_with_Bowtie.sh**: Shell script for aligning sequences using Bowtie2.
- **Trimming.sh**: Shell script for trimming adapter sequences and low-quality bases.
- **GO_condition.csv**: Results of GO Term enrichment analysis for conditions.
- **GO_time.csv**: Results of GO Term enrichment analysis for time points.
- **Quality_mapping_scores.xlsx**: Quality scores from mapping.
- **barplot_condition.pdf**: Bar plot visualization for conditions.
- **barplot_time.pdf**: Bar plot visualization for time points.
- **ego_condition.pdf**: Enriched GO Terms for conditions.
- **ego_time.pdf**: Enriched GO Terms for time points.

## How to Run the Analysis
1. **Clone the Repository**
   ```
   git clone https://github.com/yourusername/Transcriptomics_Analysis_COVID19.git
   cd Transcriptomics_Analysis_COVID19
   ```

2. **Quality Control**
   Run FastQC on your raw sequence data:
   ```sh
   fastqc raw_data/*.fastq
   ```

3. **Trimming**
   Run the trimming script:
   ```sh
   bash Trimming.sh
   ```

4. **Mapping**
   Align the trimmed sequences to the human reference genome using Bowtie2:
   ```sh
   bash Mapping_with_Bowtie.sh
   ```

5. **Differential Expression Analysis**
   Perform differential expression analysis using DESeq2:
   ```R
   Rscript Differential_expression_using_DESeq2.R
   ```

6. **GO Term Enrichment Analysis**
   Conduct GO Term enrichment analysis:
   ```R
   Rscript Gene_ontology.R
   ```

## Conclusion
This project provides a comprehensive pipeline for analyzing transcriptomic changes in human respiratory cells in response to SARS-CoV-2 infection. The integration of differential gene expression and GO Term enrichment analyses offers valuable insights into the molecular mechanisms underlying COVID-19.

## Contact
For any questions or contributions, please contact Simranjit Kang at simivk1991@gmail.com.

---

This repository provides all the necessary scripts and data files to reproduce the analysis and results presented in the study. Thank you for your interest in this project!
