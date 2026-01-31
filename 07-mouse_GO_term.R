setwd("/Users/jiaqili/Documents/Course/RNA-seq/2_Results/DESeq2")
library(grid)
library(gridExtra)
library(pathviewr) 
library(AnnotationDbi)
library(clusterProfiler)
library(topGO)
library(ggplot2)
install.packages("BiocManager")
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)

sample_name <- "DKO_vs_WT_in_infection(interative_design)"

df <- read.csv("DE_interaction_DKO_vs_WT_in_infection_with_gene_symbol.csv", sep = ",", header = T)

rownames(df) <- df$ensembl_gene_id
df <- df[order(df$padj), ]

# Define upregulated and downregulated genes based on pvalue value

genes_up <- which(df$padj < 0.05 & df$log2FoldChange > 0)
genes_down <- which(df$padj < 0.05 & df$log2FoldChange < 0)

all_genes_names <- rownames(df)

genes_up <- rownames(df)[genes_up]
genes_down <- rownames(df)[genes_down]

genelist_up <- factor(as.integer(all_genes_names %in% genes_up))
names(genelist_up) <- all_genes_names

genelist_down <- factor(as.integer(all_genes_names %in% genes_down))
names(genelist_down) <- all_genes_names

# Run topGO enrichment analysis for up- and down-regulated genes

allGO2genes <- annFUN.org(whichOnto = "ALL",
                          feasibleGenes = NULL,
                          mapping = "org.Mm.eg.db",  
                          ID = "ensembl")

GOdata_up_bp <- new("topGOdata",
                    ontology = "BP",
                    allGenes = genelist_up,
                    annot = annFUN.GO2genes,
                    GO2genes = allGO2genes,
                    nodeSize = 10)

GOdata_up_mf <- new("topGOdata",
                    ontology = "MF",
                    allGenes = genelist_up,
                    annot = annFUN.GO2genes,
                    GO2genes = allGO2genes,
                    nodeSize = 10)

GOdata_up_cc <- new("topGOdata",
                    ontology = "CC",
                    allGenes = genelist_up,
                    annot = annFUN.GO2genes, 
                    GO2genes = allGO2genes,
                    nodeSize = 10)

GOdata_down_bp <- new("topGOdata",
                      ontology = "BP",
                      allGenes = genelist_down,
                      annot = annFUN.GO2genes,
                      GO2genes = allGO2genes,
                      nodeSize = 10)

GOdata_down_mf <- new("topGOdata",
                      ontology = "MF",
                      allGenes = genelist_down,
                      annot = annFUN.GO2genes,
                      GO2genes = allGO2genes,
                      nodeSize = 10)

GOdata_down_cc <- new("topGOdata",
                      ontology = "CC",
                      allGenes = genelist_down,
                      annot = annFUN.GO2genes,
                      GO2genes = allGO2genes,
                      nodeSize = 10)


resultFis_up_bp <- runTest(GOdata_up_bp, statistic = "fisher")
resultFis_up_mf <- runTest(GOdata_up_mf, statistic = "fisher")
resultFis_up_cc <- runTest(GOdata_up_cc, statistic = "fisher")
resultFis_down_bp <- runTest(GOdata_down_bp, statistic = "fisher")
resultFis_down_mf <- runTest(GOdata_down_mf, statistic = "fisher")
resultFis_down_cc <- runTest(GOdata_down_cc, statistic = "fisher")


# Extract GO terms and p-values into a data frame
GO_terms_down_bp <- GenTable(GOdata_down_bp, pvalue = resultFis_down_bp, topNodes = length(score(resultFis_down_bp)))
GO_terms_up_bp <- GenTable(GOdata_up_bp, pvalue = resultFis_up_bp, topNodes = length(score(resultFis_up_bp)))
GO_terms_down_cc <- GenTable(GOdata_down_cc, pvalue = resultFis_down_cc, topNodes = length(score(resultFis_down_cc)))
GO_terms_up_cc <- GenTable(GOdata_up_cc, pvalue = resultFis_up_cc, topNodes = length(score(resultFis_up_cc)))

# Save the results to a CSV file
write.csv(GO_terms_down_bp, "GO_terms_down_bp.csv", row.names = FALSE)
write.csv(GO_terms_up_bp, "GO_terms_up_bp.csv", row.names = FALSE)
write.csv(GO_terms_down_cc, "GO_terms_down_cc.csv", row.names = FALSE)
write.csv(GO_terms_up_cc, "GO_terms_up_cc.csv", row.names = FALSE)

# # Extract the top 20 enriched GO terms from topGO results
parse_tables <- function(GO_data, statistics)
{
  goEnrichment <- GenTable(GO_data, weightFisher = statistics, topNodes = 20)
  sub("< ", "", goEnrichment$weightFisher)
  goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
  goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
  goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep = ", ")
  goEnrichment$Term <- factor(goEnrichment$Term, levels = rev(goEnrichment$Term))
  goEnrichment$weightFisher <- as.numeric(sub("< ", "", goEnrichment$weightFisher))  
  goEnrichment
}

#Process topGO enrichment results for upregulated and downregulated genes
# by Biological Process (BP), Molecular Function (MF), and Cellular Component (CC).
GOres_up_bp <- parse_tables(GOdata_up_bp, resultFis_up_bp)
GOres_up_mf <- parse_tables(GOdata_up_mf, resultFis_up_mf)
GOres_up_cc <- parse_tables(GOdata_up_cc, resultFis_up_cc)

GOres_down_bp <- parse_tables(GOdata_down_bp, resultFis_down_bp)
GOres_down_mf <- parse_tables(GOdata_down_mf, resultFis_down_mf)
GOres_down_cc <- parse_tables(GOdata_down_cc, resultFis_down_cc)

# Plot GO enrichment results
plot_GO <- function(GO_data, Ontology, Regulation, use_color) {
  GO_data$log_weightFisher <- (- log10(as.numeric(GO_data$weightFisher)))
  ggplot(GO_data, 
         aes(x = GO_data$log_weightFisher,
             y = GO_data$Term)) +
    geom_segment(aes(x = 0,
                     xend = GO_data$log_weightFisher,
                     y = GO_data$Term,
                     yend = GO_data$Term),
                 colour = use_color)  +
    geom_point(aes(size = GO_data$Significant),
               colour = use_color) +
    scale_size_area(name = "Gene counts") +
    xlab("Enrichment (- log10 Pvalue)") +
    ylab(Ontology) +
    ggtitle(Regulation) +
    scale_x_continuous() +
    theme_bw() +
    theme(
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_text(color = "black"))
}

plot_up_BP <- plot_GO(GOres_up_bp, "Biological Process", "TopGO Up (fisher's exact test)", "orange")
plot_up_MF <- plot_GO(GOres_up_mf, "Molecular Function", "TopGO Up (fisher's exact test)", "orange")
plot_up_CC <- plot_GO(GOres_up_cc, "Cellular Component", "TopGO Up (fisher's exact test)", "orange")

plot_down_BP <- plot_GO(GOres_down_bp, "Biological Process", "TopGO Down (fisher's exact test)", "#33CCFF")
plot_down_MF <- plot_GO(GOres_down_mf, "Molecular Function", "TopGO Down (fisher's exact test)", "#33CCFF")
plot_down_CC <- plot_GO(GOres_down_cc, "Cellular Component", "TopGO Down (fisher's exact test)", "#33CCFF")

pdf(paste(sample_name, "Biological_Proccess_TopGO_Up_fisher_pvalue.pdf", sep = "_"))
plot_up_BP
dev.off()

pdf(paste(sample_name, "Molecular_Function_TopGO_Up_fisher_pvalue.pdf", sep = "_"))
plot_up_MF
dev.off()

pdf(paste(sample_name, "Cellular_Component_TopGO_Up_fisher_pvalue.pdf", sep = "_"))
plot_up_CC
dev.off()

pdf(paste(sample_name, "Biological_Proccess_TopGO_down_fisher_pvalue.pdf", sep = "_"))
plot_down_BP
dev.off()

pdf(paste(sample_name, "Molecular_Function_TopGO_down_fisher_pvalue.pdf", sep = "_"))
plot_down_MF
dev.off()

pdf(paste(sample_name, "Cellular_Component_TopGO_down_fisher_pvalue.pdf", sep = "_"))
plot_down_CC
dev.off()

