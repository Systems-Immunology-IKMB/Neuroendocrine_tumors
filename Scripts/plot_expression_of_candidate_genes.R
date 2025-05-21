## Setup ----------------------------------------
### Bioconductor and CRAN libraries used

library(tidyverse)
library(reshape2)
library(DESeq2)
library(tximport)
library(rstatix)
library(ggpubr)
library(patchwork)

setwd("C:/Projects/Neuroendocrine_tumors")


## Loading data - organoid samples ------------------------------------------------------------

#load metadata
metadata <- read.csv("info_files/Sequenced_sample_info_organoid.csv")
metadata[metadata$Part == "Normal", ]$Group <- "control" ##combine organoids from healthy tissue independent of diagnosis of the patient for this analysis
metadata <- metadata %>%
  subset(.$RNAseq_QC == "Pass") %>%
  mutate("Part" = factor(Part, 
                         levels = c("Normal", "Tumor", "Liver metastasis")),
         "Liver.metastasis" = factor(Liver.metastasis, 
                                     levels = c("Yes", "No"))) %>%
  add_column("Part_grouped" = ifelse(.$Part == "Normal", "Normal", "Tumor") %>%
               factor(levels = c("Normal", "Tumor")))

#load count files
tx2gene <-  read.table("info_files/salmon_tx2gene.tsv")
colnames(tx2gene) <- c("Transcript_ID", "Gene_ID", "Gene_name")
annotation <- tx2gene[, 2:3] %>% unique()
rownames(annotation) <- annotation$Gene_ID

files <- list.files(path = "count_files_organoid_samples/salmon_files", full.names = T, pattern = "quant.sf", recursive = TRUE)
names(files) <- str_replace(files, "count_files_organoid_samples/salmon_files/", "") %>% 
  str_replace("_quant.sf", "") 
files <- files[metadata$Sample_ID]
txi <- tximport(files, type = "salmon", tx2gene = tx2gene[, 1:2], countsFromAbundance = "no")

all(colnames(txi$counts) == metadata$Sample_ID) #TRUE

#normalize counts
dds_counts <- DESeqDataSetFromTximport(txi, metadata, ~ Part)
dds_counts <- dds_counts[ rowSums(counts(dds_counts) == 0) < 0.8*ncol(txi$counts), ]
dds_counts <- estimateSizeFactors(dds_counts)
vst_counts <- vst(dds_counts)
vst_counts <- assay(vst_counts)

vst_counts_organoids <- as.data.frame(vst_counts)
metadata_organoids <- metadata



## Loading data - surgical samples ------------------------------------------------------------

#load metadata
metadata <- read.csv("info_files/Sequenced_sample_info_surgical.csv")
metadata <- metadata %>%
  subset(.$RNAseq_QC == "Pass") %>%
  mutate("Part" = factor(Part, 
                         levels = c("Normal", "Tumor", "Liver metastasis")),
         "Group" = factor(Group, 
                          levels = c("Insulinoma", "NF-NET")),
         "Liver.metastasis" = factor(Liver.metastasis, 
                                     levels = c("Yes", "No")))

#load count files
tx2gene <-  read.table("info_files/salmon_tx2gene.tsv")
colnames(tx2gene) <- c("Transcript_ID", "Gene_ID", "Gene_name")
annotation <- tx2gene[, 2:3] %>% unique()

files <- list.files(path = "count_files_surgical_samples/salmon_files", full.names = T, pattern = "quant.sf", recursive = TRUE)
names(files) <- str_replace(files, "count_files_surgical_samples/salmon_files/", "") %>% 
  str_replace("_quant.sf", "") 
files <- files[metadata$Sample_ID]
txi <- tximport(files, type = "salmon", tx2gene = tx2gene[, 1:2], countsFromAbundance = "no")

all(colnames(txi$counts) == metadata$Sample_ID) #TRUE

#normalize counts
dds_counts <- DESeqDataSetFromTximport(txi, metadata, ~ Part)
dds_counts <- dds_counts[ rowSums(counts(dds_counts) == 0) < 0.8*ncol(txi$counts), ]
dds_counts <- estimateSizeFactors(dds_counts)
vst_counts <- vst(dds_counts)
vst_counts <- assay(vst_counts)

vst_counts_surgical <- as.data.frame(vst_counts)
metadata_surgical <- metadata






## Plot expression of candidate genes in surgical and organoid samples as in Figure 3F and Supplementary Figure 3D --------------------

color_part <- c("Normal" = "#66C2A5", "Tumor" = "#FFAA42", "Liver metastasis" = "#DC6789")

my_genes <- c("DOCK10", "UCHL1", "CACNA2D1")
for(my_gene in my_genes){
  
  my_gene_id <- annotation %>%
    subset(.$Gene_name == my_gene)
  
  #surgical samples
  df_plot <- vst_counts_surgical[my_gene_id$Gene_ID, ] %>%
    rownames_to_column("Gene_ID") %>%
    melt(variable.name = "Sample_ID") %>%
    merge(annotation, by = "Gene_ID") %>%
    merge(metadata_surgical, by = "Sample_ID") %>%
    subset(.$Group %in% c("Insulinoma")) %>%
    subset(.$Part_grouped %in% c("Normal", "Tumor"))
  
  df_plot$Part_grouped <- factor(df_plot$Part_grouped)
  df_plot$Patient_ID <- factor(df_plot$Patient_ID)
  
  #get significance from DGE analysis
  df_res_ins <- read.csv("results/DGE_analysis/output/DEres_surgical_samples_Insulinoma_Normal_vs_Tumor.txt", sep = '\t')  %>%
    add_column("group1" = "Normal",
               "group2" = "Tumor") %>%
    rownames_to_column("Gene_ID")
  
  df_max <- df_plot %>%
    group_by(Gene_name) %>%
    dplyr::summarize("Tumor" = max(value, na.rm = TRUE) - (max(value, na.rm = TRUE) - min(value, na.rm = TRUE))*0.05) %>%
    melt(id.vars = "Gene_name", variable.name = "Group", value.name = "max")
  
  significance <- df_res_ins %>%
    subset(.$Gene_ID %in% my_gene_id$Gene_ID) %>% 
    add_significance(p.col = "padj", 
                     output.col = "p.adj_signif", 
                     cutpoints = c(0, 0.0001, 0.01, 0.05, 1),
                     symbols = c("***", "**", "*", "ns")) %>%
    merge(df_max, by.x = c("gene_name", "group2"), by.y = c("Gene_name", "Group")) %>%
    dplyr::rename("Gene_name" = "gene_name")
  
  p_surgical <- ggplot(df_plot) +
    geom_line(aes(group = Patient_ID, x = Part_grouped, y = value), position = position_dodge(width = 0.5), color = "#a9a9a9", lwd = 0.3) +
    geom_boxplot(outlier.shape = NA, aes(x = Part_grouped, y = value, color = Part_grouped, fill = Part_grouped), alpha = 0.5, lwd = 0.5) +
    geom_point(position = position_dodge(width = 0.5), aes(x = Part_grouped, y = value, color = Part_grouped, group = Patient_ID), size = 1) +
    scale_color_manual(values = color_part, breaks = c("Normal", "Tumor"), labels = c("Normal", "Insulinoma")) +
    scale_fill_manual(values = color_part, breaks = c("Normal", "Tumor"), labels = c("Normal", "Insulinoma")) +
    stat_pvalue_manual(data = significance, y.position = "max", label = "p.adj_signif", tip.length = 0, bracket.shorten = 0.5) +
    labs(x = "Insulinoma", y = "Variance stabilized counts", color = "Tissue", fill = "Tissue", title = paste0(my_gene, "\n\nSurgical samples")) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = "white"),
          axis.text.x = element_blank(),
          axis.title = element_text(size = 10),
          title = element_text(size = 10),
          legend.position = "bottom",
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank())
  
  
  
  #organoid samples
  df_plot <- vst_counts_organoids[my_gene_id$Gene_ID, ] %>%
    rownames_to_column("Gene_ID") %>%
    melt(variable.name = "Sample_ID") %>%
    merge(annotation, by = "Gene_ID") %>%
    merge(metadata_organoids, by = "Sample_ID") %>%
    subset(.$Group %in% c("control", "insulinoma")) 
  
  df_plot$Part_grouped <- factor(df_plot$Part_grouped)
  
  #get significance from DGE analysis
  df_res_ins <- read.csv("results/DGE_analysis/output/DEres_organoid_samples_Insulinoma_Normal_vs_Tumor.txt", sep = '\t')  %>%
    add_column("group1" = "Normal",
               "group2" = "Tumor") %>%
    rownames_to_column("Gene_ID")
  
  df_max <- df_plot %>%
    group_by(Gene_name) %>%
    dplyr::summarize("Tumor" = max(value, na.rm = TRUE) - (max(value, na.rm = TRUE) - min(value, na.rm = TRUE))*0.05) %>%
    melt(id.vars = "Gene_name", variable.name = "Group", value.name = "max")
  
  significance <- df_res_ins %>%
    subset(.$Gene_ID %in% my_gene_id$Gene_ID) %>% 
    add_significance(p.col = "padj", 
                     output.col = "p.adj_signif", 
                     cutpoints = c(0, 0.0001, 0.01, 0.05, 1),
                     symbols = c("***", "**", "*", "ns")) %>%
    merge(df_max, by.x = c("gene_name", "group2"), by.y = c("Gene_name", "Group")) %>%
    dplyr::rename("Gene_name" = "gene_name")
  
  p_organoid <- ggplot(df_plot) +
    geom_boxplot(outlier.shape = NA, aes(x = Part_grouped, y = value, color = Part_grouped, fill = Part_grouped), alpha = 0.5, lwd = 0.5) +
    geom_point(position = position_dodge(width = 0.5), aes(x = Part_grouped, y = value, color = Part_grouped), size = 1) +
    scale_color_manual(values = color_part, breaks = c("Normal", "Tumor"), labels = c("Normal", "Insulinoma")) +
    scale_fill_manual(values = color_part, breaks = c("Normal", "Tumor"), labels = c("Normal", "Insulinoma")) +
    stat_pvalue_manual(data = significance, y.position = "max", label = "p.adj_signif", tip.length = 0, bracket.shorten = 0.5) +
    labs(x = "Insulinoma", y = "Variance stabilized counts", color = "Tissue", fill = "Tissue", title = "\nPrimary cell\nculture samples") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = "white"),
          axis.text.x = element_blank(),
          axis.title = element_text(size = 10),
          title = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank())
  
  
  #plot expression for both surgical and organoid samples together
  my_layout <- "
AB
CC
"
  p_combined <- p_surgical + theme(legend.position = "none") +
    p_organoid + theme(legend.position = "none") +
    get_legend(p_surgical) + 
    plot_layout(heights = c(8, 2), design = my_layout)

  pdf(paste0("results/identification_of_candidate_genes/output/", my_gene, "_expression_boxplot.pdf"), width = 4, height = 4)
  plot(p_combined)
  dev.off()
  
}





