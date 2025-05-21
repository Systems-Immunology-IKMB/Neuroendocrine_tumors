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


## Loading data from surgical samples ------------------------------------------------------------

#load metadata
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

#load count data
tx2gene <-  read.table("info_files/salmon_tx2gene.tsv")
colnames(tx2gene) <- c("Transcript_ID", "Gene_ID", "Gene_name")
annotation <- tx2gene[, 2:3] %>% unique()
rownames(annotation) <- annotation$Gene_ID

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
vst_counts <- vst(dds_counts)
vst_counts <- assay(vst_counts)
vst_counts <- as.data.frame(vst_counts)

color_group <- c("NF-NET" = "#CDCDCD", "Insulinoma" = "#606060")
color_part <- c("Normal" = "#66C2A5", "Tumor" = "#FFAA42", "Liver metastasis" = "#DC6789")
color_metastasis <- c("Yes" = "#D6604D", "No" = "#4393C3")
color_grade <- c("G1" = "#9ECAE1", "G2" = "#4292C6", "G3" = "#08306B")




## Prepare insulin secretion subplot ------------------------

#import genes involved in insulin secretion (from ORA analysis)
table_rgd_pathways <- read.table("info_files/RGD_gene2pathway.txt", header = TRUE, sep = '\t', fill = TRUE, quote = "") %>%
  merge(annotation, by.x = "OBJECT_SYMBOL", by.y = "Gene_name")

df_rgd <- table_rgd_pathways[, c("Gene_ID", "OBJECT_SYMBOL", "TERM_ACC_ID", "TERM_NAME")]
colnames(df_rgd) <- c("Gene_ID", "Gene_name", "Term_ID", "Term_name")

insulin_secretion_genes <- df_rgd %>%
  subset(.$Term_ID == "PW:0000674")

#prepare data for plot
df_plot <- vst_counts[insulin_secretion_genes$Gene_ID, ] %>%
  rownames_to_column("Gene_ID") %>%
  melt(variable.name = "Sample_ID") %>%
  merge(annotation, by = "Gene_ID") %>%
  merge(metadata, by = "Sample_ID") %>%
  subset(.$Group %in% c("Insulinoma")) %>%
  subset(.$Part_grouped %in% c("Normal", "Tumor"))

df_plot$Patient_ID <- factor(df_plot$Patient_ID)
df_plot$Part_grouped <- factor(df_plot$Part_grouped)
df_plot$Gene_name <- factor(df_plot$Gene_name)

#get significance from DGE analysis
df_res_ins <- read.csv("results/DGE_analysis/output/DEres_surgical_samples_Insulinoma_Normal_vs_Tumor.txt", sep = '\t') %>%
  add_column("group1" = "Normal",
             "group2" = "Tumor") %>%
  rownames_to_column("Gene_ID")

df_max <- df_plot %>%
  group_by(Gene_name) %>%
  dplyr::summarize("Tumor" = max(value, na.rm = TRUE) - (max(value, na.rm = TRUE) - min(value, na.rm = TRUE))*0.05) %>%
  melt(id.vars = "Gene_name", variable.name = "Group", value.name = "max")

significance <- df_res_ins %>%
  subset(.$Gene_ID %in% insulin_secretion_genes$Gene_ID) %>% 
  add_significance(p.col = "padj", 
                   output.col = "p.adj_signif", 
                   cutpoints = c(0, 0.0001, 0.01, 0.05, 1),
                   symbols = c("***", "**", "*", "ns")) %>%
  merge(df_max, by.x = c("gene_name", "group2"), by.y = c("Gene_name", "Group")) %>%
  dplyr::rename("Gene_name" = "gene_name")

#prepare ggplot object
p_secretion <- ggplot(df_plot) +
  geom_line(aes(group = Patient_ID, x = Part_grouped, y = value), position = position_dodge(width = 0.2), color = "#a9a9a9", lwd = 0.3) +
  geom_boxplot(outlier.shape = NA, aes(x = Part_grouped, y = value, color = Part_grouped, fill = Part_grouped), alpha = 0.5, lwd = 0.5) +
  geom_point(position = position_dodge(width = 0.2), aes(x = Part_grouped, y = value, color = Part_grouped, group = Patient_ID), size = 1) +
  facet_wrap(. ~ Gene_name, scales = "free_y", nrow = 3) +
  scale_color_manual(values = color_part, breaks = c("Normal", "Tumor"), labels = c("Normal", "Insulinoma")) +
  scale_fill_manual(values = color_part, breaks = c("Normal", "Tumor"), labels = c("Normal", "Insulinoma")) +
  stat_pvalue_manual(data = significance, y.position = "max", label = "p.adj_signif", tip.length = 0, bracket.shorten = 0.5) +
  labs(y = "Variance stabilized counts", color = "Tissue", fill = "Tissue") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())




## Prepare insulin synthesis subplot ------------------------

#import genes involved in insulin synthesis
table_insulin_pathways <- read.table("info_files/Insulin_pathway_genes.txt", sep = "\t", header = TRUE) %>%
  merge(annotation, by = "Gene_name", all.x = TRUE)

synthesis_genes <- table_insulin_pathways %>%
  subset(.$Pathway == "Synthesis pathway")

#prepare data for plot
df_plot <- vst_counts[synthesis_genes$Gene_ID, ] %>%
  rownames_to_column("Gene_ID") %>%
  melt(variable.name = "Sample_ID") %>%
  merge(annotation, by = "Gene_ID") %>%
  merge(metadata, by = "Sample_ID") %>%
  subset(.$Group %in% c("Insulinoma")) %>%
  subset(.$Part_grouped %in% c("Normal", "Tumor"))

df_plot$Patient_ID <- factor(df_plot$Patient_ID)
df_plot$Part_grouped <- factor(df_plot$Part_grouped)
df_plot$Gene_name <- factor(df_plot$Gene_name)

#get significance from DGE analysis
df_res_ins <- read.csv("results/DGE_analysis/output/DEres_surgical_samples_Insulinoma_Normal_vs_Tumor.txt", sep = '\t') %>%
  add_column("group1" = "Normal",
             "group2" = "Tumor") %>%
  rownames_to_column("Gene_ID")

df_max <- df_plot %>%
  group_by(Gene_name) %>%
  dplyr::summarize("Tumor" = max(value, na.rm = TRUE) - (max(value, na.rm = TRUE) - min(value, na.rm = TRUE))*0.05) %>%
  melt(id.vars = "Gene_name", variable.name = "Group", value.name = "max")

significance <- df_res_ins %>%
  subset(.$Gene_ID %in% synthesis_genes$Gene_ID) %>% 
  add_significance(p.col = "padj", 
                   output.col = "p.adj_signif", 
                   cutpoints = c(0, 0.0001, 0.01, 0.05, 1),
                   symbols = c("***", "**", "*", "ns")) %>%
  merge(df_max, by.x = c("gene_name", "group2"), by.y = c("Gene_name", "Group")) %>%
  dplyr::rename("Gene_name" = "gene_name")

#prepare ggplot object
p_synthesis <- ggplot(df_plot) +
  geom_line(aes(group = Patient_ID, x = Part_grouped, y = value), position = position_dodge(width = 0.2), color = "#a9a9a9", lwd = 0.3) +
  geom_boxplot(outlier.shape = NA, aes(x = Part_grouped, y = value, color = Part_grouped, fill = Part_grouped), alpha = 0.5, lwd = 0.5) +
  geom_point(position = position_dodge(width = 0.2), aes(x = Part_grouped, y = value, color = Part_grouped, group = Patient_ID), size = 1) +
  facet_wrap(. ~ Gene_name, scales = "free_y", nrow = 3) +
  scale_color_manual(values = color_part, breaks = c("Normal", "Tumor"), labels = c("Normal", "Insulinoma")) +
  scale_fill_manual(values = color_part, breaks = c("Normal", "Tumor"), labels = c("Normal", "Insulinoma")) +
  stat_pvalue_manual(data = significance, y.position = "max", label = "p.adj_signif", tip.length = 0, bracket.shorten = 0.5) +
  labs(y = "Variance stabilized counts", color = "Tissue", fill = "Tissue") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())




## Prepare paracrine signalling subplot ------------------------

#import genes involved in paracrine signalling
table_insulin_pathways <- read.table("info_files/Insulin_pathway_genes.txt", sep = "\t", header = TRUE) %>%
  merge(annotation, by = "Gene_name", all.x = TRUE)

paracrine_genes <- table_insulin_pathways %>%
  subset(.$Pathway == "Paracrine signals")

#prepare data for plot
df_plot <- vst_counts[paracrine_genes$Gene_ID, ] %>%
  rownames_to_column("Gene_ID") %>%
  melt(variable.name = "Sample_ID") %>%
  merge(annotation, by = "Gene_ID") %>%
  merge(metadata, by = "Sample_ID") %>%
  subset(.$Group %in% c("Insulinoma")) %>%
  subset(.$Part_grouped %in% c("Normal", "Tumor"))

df_plot$Patient_ID <- factor(df_plot$Patient_ID)
df_plot$Part_grouped <- factor(df_plot$Part_grouped)
df_plot$Gene_name <- factor(df_plot$Gene_name)

#get significance from DGE analysis
df_res_ins <- read.csv("results/DGE_analysis/output/DEres_surgical_samples_Insulinoma_Normal_vs_Tumor.txt", sep = '\t') %>%
  add_column("group1" = "Normal",
             "group2" = "Tumor") %>%
  rownames_to_column("Gene_ID")

df_max <- df_plot %>%
  group_by(Gene_name) %>%
  dplyr::summarize("Tumor" = max(value, na.rm = TRUE) - (max(value, na.rm = TRUE) - min(value, na.rm = TRUE))*0.05) %>%
  melt(id.vars = "Gene_name", variable.name = "Group", value.name = "max")

significance <- df_res_ins %>%
  subset(.$Gene_ID %in% paracrine_genes$Gene_ID) %>% 
  add_significance(p.col = "padj", 
                   output.col = "p.adj_signif", 
                   cutpoints = c(0, 0.0001, 0.01, 0.05, 1),
                   symbols = c("***", "**", "*", "ns")) %>%
  merge(df_max, by.x = c("gene_name", "group2"), by.y = c("Gene_name", "Group")) %>%
  dplyr::rename("Gene_name" = "gene_name")

#prepare ggplot object
p_paracrine <- ggplot(df_plot) +
  geom_line(aes(group = Patient_ID, x = Part_grouped, y = value), position = position_dodge(width = 0.2), color = "#a9a9a9", lwd = 0.3) +
  geom_boxplot(outlier.shape = NA, aes(x = Part_grouped, y = value, color = Part_grouped, fill = Part_grouped), alpha = 0.5, lwd = 0.5) +
  geom_point(position = position_dodge(width = 0.2), aes(x = Part_grouped, y = value, color = Part_grouped, group = Patient_ID), size = 1) +
  facet_wrap(. ~ Gene_name, scales = "free_y", nrow = 3) +
  scale_color_manual(values = color_part, breaks = c("Normal", "Tumor"), labels = c("Normal", "Insulinoma")) +
  scale_fill_manual(values = color_part, breaks = c("Normal", "Tumor"), labels = c("Normal", "Insulinoma")) +
  stat_pvalue_manual(data = significance, y.position = "max", label = "p.adj_signif", tip.length = 0, bracket.shorten = 0.5) +
  labs(y = "Variance stabilized counts", color = "Tissue", fill = "Tissue") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())




## Expression boxplots for genes in different insulin pathways as in Figure 2D ----------------

p_combined <- p_synthesis + labs(title = "Synthesis\npathway") + theme(legend.position = "none") +
  p_secretion + labs(title = "Secretion\npathway") + theme(axis.title.y = element_blank()) +
  p_paracrine + labs(title = "Paracrine\npathway") + theme(legend.position = "none", axis.title.y = element_blank()) +
  plot_layout(nrow = 1, width = c(7, 6, 0.8))

pdf("results/DGE_analysis/output/expression_of_insulin_pathway_genes.pdf", width = 21, height = 7)
p_combined
dev.off()



