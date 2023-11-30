library(DESeq2)

# Read in count data and metadata
counts <- read.csv("counts_for_deseq2.csv", row.names=1)
metadata <- read.csv("metadata_for_deseq2.csv")

# Create DESeq dataset
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ condition + time)

# Run the DESeq pipeline
dds <- DESeq(dds)

# Get results for condition comparison
res <- results(dds, contrast=c("condition", "Control", "Treatment"))

# Export results
write.csv(as.data.frame(res), file="DESeq2_results.csv")

cat("DESeq2 analysis complete.\n")
