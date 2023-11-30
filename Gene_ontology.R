#!/usr/bin/env Rscript

library(clusterProfiler)

# Assuming 'res' contains the DESeq2 results with a column 'gene' for gene IDs
res <- read.csv("DESeq2_results.csv")

# Extract significantly differentially expressed genes
sig_genes <- res$gene[res$padj < 0.05]

# Perform GO enrichment analysis
ego <- enrichGO(gene = sig_genes,
                OrgDb = org.Hs.eg.db,
                keyType = "ENSEMBL",
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)

# Export GO enrichment analysis results
write.csv(as.data.frame(ego), file="GO_enrichment_results.csv")

cat("GO enrichment analysis complete.\n")

