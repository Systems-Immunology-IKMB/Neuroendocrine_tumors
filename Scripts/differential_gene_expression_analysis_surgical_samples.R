## Setup ----------------------------------------
### Bioconductor and CRAN libraries used

library(tidyverse)
library(reshape2)
library(DESeq2)
library(tximport)
library(ggrepel)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

setwd("C:/Projects/Neuroendocrine_tumors")


## Loading data ------------------------------------------------------------

metadata <- read.csv("info_files/Sequenced_sample_info_surgical.csv")
metadata <- metadata %>%
  subset(.$RNAseq_QC == "Pass") %>%
  mutate("Part" = factor(Part, 
                         levels = c("Normal", "Tumor", "Liver metastasis")),
         "Group" = factor(Group, 
                          levels = c("Insulinoma", "NF-NET")),
         "Liver.metastasis" = factor(Liver.metastasis, 
                                     levels = c("Yes", "No"))) %>%
  add_column("Part_grouped" = ifelse(.$Part == "Normal", "Normal", "Tumor") %>%
               factor(levels = c("Normal", "Tumor")))

tx2gene <-  read.table("info_files/salmon_tx2gene.tsv")
colnames(tx2gene) <- c("Transcript_ID", "Gene_ID", "Gene_name")
annotation <- tx2gene[, 2:3] %>% unique()
rownames(annotation) <- annotation$Gene_ID



## Differential gene expression analysis - normal vs insulinoma ----------------------

#prepare metadata
sample_info_filtered <- subset(metadata, metadata$Group == "Insulinoma" & metadata$Part_grouped %in% c("Normal", "Tumor"))

sample_info_filtered$Patient_ID <- as.factor(sample_info_filtered$Patient_ID)
sample_info_filtered$Grade <- as.factor(sample_info_filtered$Grade)
sample_info_filtered$Part_grouped <- as.factor(sample_info_filtered$Part_grouped)

#import count data
files <- paste("count_files_surgical_samples/salmon_files/", sample_info_filtered$Sample_ID, "_quant.sf", sep = '')
names(files) <- sample_info_filtered$Sample_ID
txi <- tximport(files, type = "salmon", tx2gene = tx2gene[, 1:2], countsFromAbundance = "no")
print(all(colnames(txi$counts) == sample_info_filtered$Sample_ID)) #TRUE

#run DE analysis
dds_counts <- DESeqDataSetFromTximport(txi, sample_info_filtered, ~ Patient_ID + Part_grouped)
dds_counts <- dds_counts[ rowSums(counts(dds_counts) == 0) < 0.8*ncol(txi$counts), ]
dds_counts <- estimateSizeFactors(dds_counts)
dds <- DESeq(dds_counts, betaPrior = FALSE)
res <- results(dds, independentFiltering = TRUE, alpha = 0.05)
res_sorted <- res[order(res$pvalue), ]
res_sorted$gene_name <- annotation[rownames(res_sorted), "Gene_name"]

output_directory <- paste0("results/DGE_analysis/output/")
write.table(res_sorted, file = paste0(output_directory, "DEres_surgical_samples_Insulinoma_Normal_vs_Tumor.txt"), sep = "\t", quote = FALSE)





## Volcano plot for insulinoma DEGs as in Figure 2A ------------------------------------------------------------

#load DGE analysis results and filter for DEGs
res_table_insulinoma <- read.csv("results/DGE_analysis/output/DEres_surgical_samples_Insulinoma_Normal_vs_Tumor.txt", sep = '\t')
res_table_insulinoma <- res_table_insulinoma %>%
  subset(.$padj < 0.05) %>%
  subset(.$baseMean > 10) %>%
  subset(abs(.$log2FoldChange) > 1) %>%
  rownames_to_column("Gene_ID")

df_plot <- res_table_insulinoma
df_plot$log_padj <- -1 * log10(df_plot$padj)
df_plot$color <- ifelse(df_plot$log2FoldChange > 1 & df_plot$padj < 0.05, "red", 
                        ifelse(df_plot$log2FoldChange < -1 & df_plot$padj < 0.05, "blue", 
                               ifelse(df_plot$padj < 0.05, "black", "grey")))

#add labels for top genes
label_fc <- df_plot %>% 
  .[order(.$log2FoldChange, decreasing = TRUE), ] %>%
  {c(head(., 5)$Gene_ID, tail(., 2)$Gene_ID)}

label_pval_up <- df_plot %>% 
  subset(.$log2FoldChange > 0) %>%
  .[order(.$log_padj, decreasing = TRUE), ] %>%
  {c(head(., 2)$Gene_ID)}

label_pval_down <- df_plot %>% 
  subset(.$log2FoldChange < 0) %>%
  .[order(.$log_padj, decreasing = TRUE), ] %>%
  {c(head(., 5)$Gene_ID)}

df_plot$label <- ifelse(df_plot$Gene_ID %in% c(label_fc, label_pval_up, label_pval_down), df_plot$gene_name, NA)

#volcano plot for insulinoma DEGs
pdf("results/DGE_analysis/output/volcanoplot_insulinoma_surgical_samples.pdf", width = 8, height = 5)
ggplot(data = df_plot, aes(x = log2FoldChange, y = log_padj, color = color)) + 
  geom_point(pch = 20, size = 3) +
  geom_text_repel(aes(label = label, y = log_padj), 
                  size = 9/.pt, 
                  min.segment.length = 0, 
                  max.overlaps = 5) +
  geom_vline(xintercept = 0, lty = 3) + 
  geom_hline(yintercept = 0, lty = 3) +
  geom_vline(xintercept = 1, lty = 4) + 
  geom_vline(xintercept = -1, lty = 4) +
  geom_hline(yintercept = 1.30103, lty = 4) + 
  scale_color_manual(values = c("black"='#606060',"blue"='#3399FF',"grey"='#C0C0C0',"red"='#FF6666')) +
  xlab("Expression (Log2 Fold Change)") + 
  ylab("-log10 q-value") +
  theme_bw() + 
  theme(axis.text = element_text(size = 14), 
        axis.title = element_text(size = 16), 
        panel.grid = element_blank(),
        legend.position = "none")
dev.off()

df_plot %>% 
  add_column("direction" = ifelse(.$log2FoldChange > 0, "up", "down")) %>%
  {table(.$direction)}
#down   up 
#1192 2141





## Heatmap for insulinoma DEGs in normal, insulinoma and NF-NET surgical samples as in Figure 2B -------------------

#load count data
files <- list.files(path = "count_files_surgical_samples/salmon_files", full.names = T, pattern = "quant.sf", recursive = TRUE)
names(files) <- str_replace(files, "count_files_surgical_samples/salmon_files/", "") %>% 
  str_replace("_quant.sf", "") 
files <- files[metadata$Sample_ID]
txi <- tximport(files, type = "salmon", tx2gene = tx2gene[, 1:2], countsFromAbundance = "no")
all(colnames(txi$counts) == metadata$Sample_ID) #TRUE

#normalize expression counts
dds_counts <- DESeqDataSetFromTximport(txi, metadata, ~ Part)
dds_counts <- dds_counts[ rowSums(counts(dds_counts) == 0) < 0.8*ncol(txi$counts), ]
dds_counts <- estimateSizeFactors(dds_counts)
normalized_counts <- counts(dds_counts, normalized = TRUE) %>%
  as.data.frame()

#load and identify DEGs to plot
res_table_insulinoma <- read.csv("results/DGE_analysis/output/DEres_surgical_samples_Insulinoma_Normal_vs_Tumor.txt", sep = '\t')
res_table_insulinoma <- res_table_insulinoma %>%
  subset(.$padj < 0.05) %>%
  subset(.$baseMean > 10) %>%
  subset(abs(.$log2FoldChange) > 1) %>%
  rownames_to_column("Gene_ID")
gene_ids <- res_table_insulinoma$Gene_ID
degs <- annotation[annotation$Gene_ID %in% gene_ids, "Gene_ID"] %>%
  subset(. %in% rownames(normalized_counts))

#prepare metadata and counts for heatmap
meta_heat <- metadata %>%
  .[order(.$Group, .$Part_grouped), ]

norm_counts_temp <- normalized_counts %>%
  .[degs, meta_heat$Sample_ID] %>%
  as.data.frame()

matDataHeat <- t(scale(t(norm_counts_temp)))

#add column annotation
color_group <- c("NF-NET" = "#CDCDCD", "Insulinoma" = "#606060")
color_part <- c("Normal" = "#66C2A5", "Tumor" = "#FFAA42", "Liver metastasis" = "#DC6789")
color_metastasis <- c("Yes" = "#D6604D", "No" = "#4393C3")
color_grade <- c("G1" = "#9ECAE1", "G2" = "#4292C6", "G3" = "#08306B")

col_anno <- meta_heat %>%
  .[match(colnames(norm_counts_temp), meta_heat$Sample_ID), c("Group", "Part", "Grade", "Liver.metastasis")]

col_ha <- HeatmapAnnotation(df = col_anno,
                            col = list(Group = color_group, 
                                       Part = color_part,
                                       Grade = color_grade,
                                       Liver.metastasis = color_metastasis), 
                            annotation_label = c("Group", "Tissue", "Grade", "Liver metastasis"),
                            border = TRUE)

col_split <- split(meta_heat, ~ meta_heat$Part_grouped + meta_heat$Group) %>%
  map_int(function(x){nrow(x)}) %>%
  {factor(rep(names(.), .), levels = names(.))} 

pdf(paste0("results/DGE_analysis/output/heatmap_insulinoma_DEGs_surgical_samples.pdf"), width = 6, height = 8)
p <- Heatmap(matDataHeat, 
             name = "Scaled\nnormalized\nexpression",
             col=colorRamp2(seq(from = -2.1, to = 2.1, by = 0.42) ,rev(brewer.pal(11, "RdYlBu"))),
             cluster_columns = TRUE,
             cluster_rows = TRUE,
             top_annotation = col_ha, 
             row_names_gp = gpar(fontsize = 7),
             column_names_gp = gpar(fontsize = 10),
             row_gap = unit(1, "mm"), 
             border = TRUE, 
             show_row_names = FALSE,
             show_row_dend = FALSE,
             show_column_names = FALSE,
             row_title_gp = gpar(fontsize = 10),
             row_title_rot = 0,
             cluster_row_slices = FALSE,
             column_split = col_split,
             column_title = NULL,
             show_column_dend = TRUE,
             cluster_column_slices = FALSE)
draw(p, merge_legend = TRUE)
dev.off()






