setwd("/Users/jiaqili/Documents/Course/RNA-seq/2_Results/DESeq2")

library(tidyverse) 
library(RColorBrewer) 
library(ggrepel) 

df <- read.csv("DE_interaction_DKO_vs_WT_in_infection_with_gene_symbol.csv", sep =",", header = T)

# Replace 0 values with the smallest non-zero padj value.
min_padj <- min(df$padj[df$padj > 0], na.rm = TRUE)
df$padj[df$padj == 0] <- min_padj

# =====================================
#definition of up-regulated and down-regulated genes
df$diffexpressed <- "NO"
df$diffexpressed[df$log2FoldChange > 0.415037 & df$padj < 0.05] <- "UP"
df$diffexpressed[df$log2FoldChange < -0.415037 & df$padj < 0.05] <- "DOWN" 

#define fow which genes genenames are displayed

# Define which genes to label: top 20 by lowest padj + top 5 by absolute log2FC
top_padj_genes <- head(df[order(df$padj), "gene_label"], 20)


top_logfc_genes <- head(df[order(-abs(df$log2FoldChange)), "gene_label"], 5)
highlight_genes <- unique(c(top_padj_genes, top_logfc_genes))

df$delabel <- ifelse(df$gene_label %in% highlight_genes, df$gene_label, NA)


theme_set(theme_classic(base_size = 20) +  #general layout (font size etc.)
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.0), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.0), color = 'black'),
              plot.title = element_text(hjust = 0.5)
              ))
ggplot(data = df, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-0.415037, 0.415037), col = "gray", linetype = 'dashed') +  #threshold lines
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + #threshold lines
  geom_point(size = 2) +
  scale_color_manual(values = c("#33CCFF", "grey", "#FF9900"), #color of up and downregulated
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  coord_cartesian(ylim = c(0, 550), xlim = c(-8, 8)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = '', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"padj")) + 
  scale_x_continuous(breaks = seq(-8, 8, 1)) + # to customise the breaks in the x axis
  ggtitle('DKO_vs_WT_in_infection(interactive_design)') + # Plot title 
  geom_text_repel(max.overlaps = Inf) # To show all labels 
geom_point()

