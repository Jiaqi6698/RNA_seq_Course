
setwd("/Users/jiaqili/Documents/Course/RNA-seq/2_Results")
OUT_DIR <- "DESeq2"
dir.create(OUT_DIR)
suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(pheatmap)
})

## Define the correspondence between samples and groups.
sample_mapping <- data.frame(
  SRR_ID = c("SRR7821921", "SRR7821922", "SRR7821918", "SRR7821919", "SRR7821920",
             "SRR7821937", "SRR7821938", "SRR7821939",
             "SRR7821923", "SRR7821924", "SRR7821925", "SRR7821927",
             "SRR7821940", "SRR7821941", "SRR7821942"),
  Sample_Name = c("Lung_WT_Case", "Lung_WT_Case", "Lung_WT_Case", "Lung_WT_Case", "Lung_WT_Case",
                  "Lung_WT_Control", "Lung_WT_Control", "Lung_WT_Control",
                  "Lung_DKO_Case", "Lung_DKO_Case", "Lung_DKO_Case", "Lung_DKO_Case",
                  "Lung_DKO_Control", "Lung_DKO_Control", "Lung_DKO_Control"),
  stringsAsFactors = FALSE
)

##Automatically add serial numbers to duplicate Sample_Name to generate unique names New_Name
sample_mapping$New_Name <- ave(sample_mapping$Sample_Name,
                               sample_mapping$Sample_Name,
                               FUN = function(x) {
                                 paste0(x, "_", seq_along(x))
                               })

## Read count matrix.
counts_data <- read.delim("counts_matrix.txt",
                          header = TRUE,
                          row.names = 1,
                          comment.char = "#",
                          check.names = FALSE)

## The counts of the sample start from column 7, so select 7:ncol
counts_matrix <- counts_data[, 7:ncol(counts_data)]

## keep the old column names for parsing SRR numbers
old_colnames <- colnames(counts_matrix)

## Extract SRR ID from column names
srr_ids <- sub(".*(SRR[0-9]+).*", "\\1", old_colnames)

## Rearrange according to the position of SRR ID in sample_mapping to align sample information.
idx <- match(srr_ids, sample_mapping$SRR_ID)
if (any(is.na(idx))) {
  stop("One or more SRR IDs from counts_matrix do not exist in sample_mapping")
}

## Generate a unique new sample name (e.g., Lung_WT_Case_1) and assign it to counts_matrix
new_names <- sample_mapping$New_Name[idx]
colnames(counts_matrix) <- new_names

## Construct the sample information table according to the rearranged idx.
## Determine group, genotype, and condition based on Sample_Name in sample_mapping.
sample_info <- data.frame(
  sample    = new_names,
  group     = sample_mapping$Sample_Name[idx],
  genotype  = ifelse(grepl("WT", sample_mapping$Sample_Name[idx]), "WT", "DKO"),
  condition = ifelse(grepl("Case", sample_mapping$Sample_Name[idx]), "Case", "Control"),
  row.names = new_names,
  stringsAsFactors = FALSE
)

## Combining genotype and condition facilitates subsequent visualization.
sample_info$genotype_condition <- factor(
  paste(sample_info$genotype, sample_info$condition, sep = "_")
)

## ---- set factor----
sample_info$genotype  <- factor(sample_info$genotype,
                                levels = c("WT", "DKO"))
sample_info$condition <- factor(sample_info$condition,
                                levels = c("Control", "Case"))

## ---- Construct DESeq2 object (interaction design)----
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData   = sample_info,
  design    = ~ genotype + condition + genotype:condition)

# Filter lowly expressed genes
dds <- dds[rowSums(counts(dds)) >= 10, ]

## ----  DESeq2 ----
dds <- DESeq(dds)

saveRDS(dds, file = file.path(OUT_DIR, "dds_full.rds"))

## ---- View coefficient name----
resultsNames(dds)


BiocManager::install("apeglm")
library("apeglm")


## ---- 10. Interaction Effect: Impact of DKO on Infection Response ----
res_interaction <- results(dds,
                           name = "genotypeDKO.conditionCase")

res_interaction <- lfcShrink(dds,
                             coef = "genotypeDKO.conditionCase",
                             res = res_interaction)

write.csv(as.data.frame(res_interaction),
          file = file.path(OUT_DIR,
                           "DE_interaction_DKO_vs_WT_in_infection.csv"))

## ---- PCA: Based on vst transformation ----
vsd <- vst(dds, blind = FALSE)
saveRDS(vsd, file = file.path(OUT_DIR, "vsd.rds"))

## ========== Plot PCA ==========
pcaData <- plotPCA(vsd,
                   intgroup = c("genotype", "condition"),
                   returnData = TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

p_pca <- ggplot(pcaData, aes(x = PC1, y = PC2,
                             color = genotype,
                             shape = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  coord_fixed() +
  theme_classic()

ggsave(filename = file.path(OUT_DIR, "PCA_plotPCA_PC1_PC2.pdf"),
       plot = p_pca, width = 6, height = 5)
ggsave(filename = file.path(OUT_DIR, "PCA_plotPCA_PC1_PC2.png"),
       plot = p_pca, width = 6, height = 5, dpi = 300)

write.csv(pcaData, file = file.path(OUT_DIR, "PCA_plotPCA_scores.csv"),
          row.names = FALSE)
