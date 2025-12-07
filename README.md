# RNAseq_pipeline
RNA-seq pipeline for DEA

###############################################################
# RNA-Seq Differential Expression and Network Analysis Pipeline
# -------------------------------------------------------------
# Author: Paraskevi Tsortanidou
# Thesis Appendix
###############################################################

# ---------------------------
# 0. Package Installation
# ---------------------------
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")

BiocManager::install(c("DESeq2", "apeglm", "vsn", "pheatmap",
                       "EnhancedVolcano", "clusterProfiler",
                       "org.Hs.eg.db", "enrichplot", "GENIE3"))

install.packages(c("readxl", "openxlsx", "ggplot2", "dplyr"))

# ---------------------------
# 1. Load Libraries
# ---------------------------
library(DESeq2)
library(apeglm)
library(vsn)
library(pheatmap)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(GENIE3)
library(readxl)
library(openxlsx)
library(ggplot2)
library(dplyr)

# ---------------------------
# 2. Load Data
# ---------------------------
counts_df  <- read_excel("merged_counts.xlsx")
coldata    <- read_excel("coldata.xlsx")

annot <- counts_df %>% select(Symbol, Name)
count_cols <- setdiff(colnames(counts_df), c("Symbol", "Name"))
count_mat <- as.matrix(counts_df[, count_cols])
rownames(count_mat) <- counts_df$Symbol
count_mat <- round(apply(count_mat, 2, as.numeric))
storage.mode(count_mat) <- "integer"

# Match metadata
coldata <- coldata[match(colnames(count_mat), coldata$sample), ]
coldata$genotype <- factor(coldata$genotype, levels = c("WT", "KO"))
coldata$treatment <- factor(coldata$treatment, levels = c("ctrl", "GHRH15", "GHRH30"))

# ---------------------------
# 3. Differential Expression (DESeq2)
# ---------------------------
dds <- DESeqDataSetFromMatrix(countData = count_mat,
                              colData = as.data.frame(coldata),
                              design = ~ genotype + treatment + genotype:treatment)

# Filter and run
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
dds <- DESeq(dds)

# ---------------------------
# 4. PCA and Sample Clustering
# ---------------------------
vsd <- vst(dds, blind = FALSE)

## PCA Plot
pcaData <- plotPCA(vsd, intgroup = c("genotype", "treatment"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = genotype, shape = treatment)) +
  geom_point(size = 4) +
  labs(
    title = "Principal Component Analysis",
    x = paste0("PC1 (", percentVar[1], "% variance)"),
    y = paste0("PC2 (", percentVar[2], "% variance)")
  ) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

## Sample Distance Heatmap
sampleDistMat <- as.matrix(dist(t(assay(vsd))))
rownames(sampleDistMat) <- colnames(vsd)
ann_col <- data.frame(Genotype = coldata$genotype,
                      Treatment = coldata$treatment)
rownames(ann_col) <- coldata$sample

pheatmap(sampleDistMat,
         annotation_col = ann_col,
         main = "Sample-to-Sample Distance (VST)",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         fontsize = 12, border_color = NA)

# ---------------------------
# 5. DE Results & Volcano Plot
# ---------------------------
res_genotype_ctrl <- results(dds, name = "genotype_KO_vs_WT", alpha = 0.05)
res_genotype_ctrl_shr <- lfcShrink(dds, coef = "genotype_KO_vs_WT", type = "apeglm")

res_df <- as.data.frame(res_genotype_ctrl_shr)
res_df$Symbol <- rownames(res_df)
res_df <- left_join(res_df, annot, by = "Symbol")

EnhancedVolcano(
  res_df,
  lab = ifelse(is.na(res_df$Name), res_df$Symbol, res_df$Name),
  x = 'log2FoldChange', y = 'padj',
  title = 'KO vs WT (control)',
  subtitle = 'Differential Expression (DESeq2, FDR < 0.05)',
  pCutoff = 0.05, FCcutoff = 1.0,
  col = c('grey70','steelblue2','gold2','firebrick')
)

# ---------------------------
# 6. Heatmap of Top 50 DE Genes
# ---------------------------
res_df <- res_df[order(res_df$padj), ]
top50 <- rownames(res_df)[1:50]
mat_top50 <- assay(vsd)[top50, ]
mat_top50 <- t(scale(t(mat_top50)))  # Z-score normalization

pheatmap(mat_top50,
         annotation_col = ann_col,
         main = "Top 50 Differentially Expressed Genes\n(KO vs WT)",
         fontsize_row = 6,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         border_color = NA)

# ---------------------------
# 7. KEGG Pathway Enrichment
# ---------------------------
dir.create("KEGG_results", showWarnings = FALSE)

run_kegg <- function(res_object, contrast_name) {
  degs <- rownames(res_object)[which(res_object$padj < 0.05 & abs(res_object$log2FoldChange) > 1)]
  if (length(degs) == 0) return(NULL)
  entrez_ids <- na.omit(mapIds(org.Hs.eg.db, keys = gsub("\\..*", "", degs),
                               keytype = "ENSEMBL", column = "ENTREZID", multiVals = "first"))
  kegg_res <- enrichKEGG(gene = entrez_ids, organism = "hsa",
                         pAdjustMethod = "BH", qvalueCutoff = 0.05)
  write.xlsx(as.data.frame(kegg_res),
             paste0("KEGG_results/KEGG_", contrast_name, ".xlsx"))
  dotplot(kegg_res, showCategory = 15) +
    ggtitle(paste("KEGG Pathway Enrichment:", contrast_name)) +
    theme_minimal(base_size = 12)
}

kegg_KO_WT <- run_kegg(res_genotype_ctrl_shr, "KO_vs_WT")

# ---------------------------
# 8. GENIE3 Gene Regulatory Network Inference
# ---------------------------
expr_matrix <- assay(vsd)
deg_genes <- rownames(res_df[res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, ])
expr_matrix_deg <- expr_matrix[deg_genes, ]

set.seed(123)
weight_matrix <- GENIE3(expr_matrix_deg)
link_list <- getLinkList(weight_matrix)
top_links <- link_list[1:500, ]

library(igraph)
g <- graph_from_data_frame(top_links, directed = TRUE)
plot(g,
     layout = layout_with_fr(g),
     vertex.size = 5, vertex.label = NA,
     edge.arrow.size = 0.3,
     main = "Top 500 Predicted Regulatory Interactions (GENIE3)")

# Export for Cytoscape
write.csv(link_list, "GENIE3_links_for_Cytoscape.csv", row.names = FALSE)
