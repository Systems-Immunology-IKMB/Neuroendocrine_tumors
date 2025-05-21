## Setup ----------------------------------------
### Bioconductor and CRAN libraries used

library(tidyverse)
library(reshape2)
library(VennDiagram)
library(dplyr)
library(ggpubr)
library(ggrepel)
library(tximport)
library(DESeq2)
library(nlme)

setwd("C:/Projects/Neuroendocrine_tumors")


## Identify genes dysregulated in insulinoma vs normal in both organoid and surgical samples -------------

#load DEGs
surgical_samples_result <- read.csv("results/DGE_analysis/output/DEres_surgical_samples_Insulinoma_Normal_vs_Tumor.txt", sep = '\t')
surgical_samples_result <- surgical_samples_result %>%
  subset(.$padj < 0.05) %>%
  subset(.$baseMean > 10) %>%
  subset(abs(.$log2FoldChange) > 1) %>%
  rownames_to_column("Gene_ID")

organoid_result <- read.csv("results/DGE_analysis/output/DEres_organoid_samples_Insulinoma_Normal_vs_Tumor.txt", sep = '\t')
organoid_result <- organoid_result %>%
  subset(.$padj < 0.05) %>%
  subset(.$baseMean > 10) %>%
  subset(abs(.$log2FoldChange) > 1) %>%
  rownames_to_column("Gene_ID")

#merge and export results
shared_degs <- inner_join(surgical_samples_result, organoid_result, by = c("Gene_ID", "gene_name"), suffix = c(".Surgical", ".Organoid"))
shared_degs$Direction_in_Surgical <- ifelse(shared_degs$log2FoldChange.Surgical > 0, "Upregulated in insulinoma", "Downregulated in insulinoma")
shared_degs$Direction_in_Organoid <- ifelse(shared_degs$log2FoldChange.Organoid > 0, "Upregulated in insulinoma", "Downregulated in insulinoma")
write.csv(shared_degs, "results/identification_of_candidate_genes/output/shared_insulinoma_DEGs_organoid_and_surgical.csv", quote = FALSE, row.names = FALSE)

#plot Venn diagram as in Figure 3A
pdf("results/identification_of_candidate_genes/output/shared_insulinoma_DEGs_vennDiagram.pdf", width = 6, height = 6)
grid.newpage()
draw.pairwise.venn(area1 = nrow(surgical_samples_result), 
                   area2 = nrow(organoid_result), 
                   cross.area = nrow(shared_degs), 
                   category = c("Surgical samples", "Cell culture\nsamples"), 
                   col = c("#024E37", "#01587A"), 
                   fill = c("#026D4E", "#077DAA"), 
                   alpha = rep(0.5, 2), 
                   cat.pos = c(0, 0), 
                   cat.dist = c(0.05, 0.05), 
                   cat.cex = rep(1.7, 2), 
                   cat.fontfamily = rep("Helvetica", 2), 
                   cex = rep(1.5, 3), 
                   fontfamily = rep("Helvetica", 3))
dev.off()


#compare LFC between organoid and surgical samples in a dotplot as in Figure 3B
df_plot <- shared_degs
label_degs_org <- df_plot %>% 
  .[order(.$log2FoldChange.Organoid, decreasing = TRUE), ] %>%
  {c(head(., 5)$Gene_ID, tail(., 10)$Gene_ID)}
label_degs_sur <- df_plot %>% 
  .[order(.$log2FoldChange.Surgical, decreasing = TRUE), ] %>%
  {c(head(., 5)$Gene_ID, tail(., 10)$Gene_ID)}
label_diff <- df_plot %>% 
  subset(.$log2FoldChange.Organoid > 0 & .$log2FoldChange.Surgical < 0 ) %>%
  .[order(.$log2FoldChange.Organoid - .$log2FoldChange.Surgical, decreasing = TRUE), ] %>%
  {head(., 3)$Gene_ID}

df_plot$label <- ifelse(df_plot$Gene_ID %in% c(label_degs_org, label_degs_sur, label_diff, "ENSG00000135905", "ENSG00000153956"), df_plot$gene_name, NA)

pdf("results/identification_of_candidate_genes/output/shared_insulinoma_DEGs_LFC_dotplot.pdf", width = 6, height = 6)
ggplot(df_plot, aes(x = log2FoldChange.Surgical, y = log2FoldChange.Organoid)) +
  geom_point(size = 2, color = "#439692", alpha = 0.5) +
  stat_cor(method = "spearman", cor.coef.name = "rho") +
  geom_text_repel(aes(label = label), size = 8/.pt, min.segment.length = 0) +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) +
  labs(x = "LFC surgical samples", y = "LFC cell culture samples") +
  theme_bw() + 
  theme(axis.text = element_text(size = 10, color = "black"), 
        panel.grid = element_blank(),
        axis.title = element_text(size = 11))
dev.off()






## Identify genes associated with insulin expression ---------------

#load metadata surgical samples
metadata <- read.csv("info_files/Sequenced_sample_info_surgical.csv")
metadata <- metadata %>%
  subset(.$RNAseq_QC == "Pass") %>%
  mutate("Part" = factor(Part, 
                         levels = c("Normal", "Tumor", "Liver metastasis")),
         "Group" = factor(Group, 
                          levels = c("Insulinoma", "NF-NET")),
         "Liver.metastasis" = factor(Liver.metastasis, 
                                     levels = c("Yes", "No")))

#load count data surgical samples
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


#normalize expression counts and subset to insulinoma samples
dds_counts <- DESeqDataSetFromTximport(txi, metadata, ~ Part)
dds_counts <- dds_counts[ rowSums(counts(dds_counts) == 0) < 0.8*ncol(txi$counts), ]
dds_counts <- estimateSizeFactors(dds_counts)
vst_counts <- vst(dds_counts)
vst_counts <- assay(vst_counts)
vst_counts <- as.data.frame(vst_counts)

metadata_insulinoma <- metadata %>%
  subset(.$Group %in% c("Insulinoma")) %>%
  subset(.$Part %in% c("Tumor", "Liver metastasis"))
vst_counts <- vst_counts[, metadata_insulinoma$Sample_ID]


#load DEGs upregulated in insulinoma and shared between organoid and surgical analysis
group <- "insulinoma"
shared_degs <- read.table("results/identification_of_candidate_genes/output/shared_insulinoma_DEGs_organoid_and_surgical.csv", header = TRUE)
shared_degs_up <- shared_degs %>%
  subset(.$Direction_in_Surgical == .$Direction_in_Organoid) %>%
  subset(.$Direction_in_Surgical == "Upregulated")
my_genes <- shared_degs_up$Gene_ID 

#Estimate association between candidate genes and INS expression
res_table <- data.frame()
for(i in 1:length(my_genes)){
  
  test_counts <- vst_counts[c("ENSG00000254647", my_genes[i]), ] %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Sample_ID") %>%
    merge(metadata)
  colnames(test_counts)[2:3] <- c("INS", "test_gene")
  
  H_1 <- lme(fixed = INS ~ test_gene, random = ~ 1|Patient_ID, data = test_counts, method = "ML")
  res_LMM <- summary(H_1)$tTable %>% as.data.frame()
  
  res_table <- rbind(res_table, data.frame("Gene_ID" = my_genes[i], 
                                               "p_value" = res_LMM$`p-value`[2],
                                               "estimate" = res_LMM$Value[2], 
                                               "t-value" = res_LMM$`t-value`[2]))
}

res_table <- res_table %>%
  merge(annotation, by = "Gene_ID") %>%
  add_column("p_adj" = p.adjust(.$p_value, method = "BH")) %>%
  .[order(.$p_value, decreasing = FALSE), ]

write.csv(res_table, "results/identification_of_candidate_genes/output/shared_DEGs_LMM_association_with_INS.csv", row.names = FALSE)





## Plot LMM results ---------------

res_table <- read.csv("results/identification_of_candidate_genes/output/shared_DEGs_LMM_association_with_INS.csv") %>%
  subset(.$estimate > 0) %>%
  add_column("log10p" = -log10(.$p_value)) 


#barplot for top 10 genes associated with INS expression as in Figure 3C
df_barplot <- res_table %>%
  .[order(.$log10p, decreasing = TRUE), ] %>%
  .[1:10, ]
df_barplot$Gene_name <- factor(df_barplot$Gene_name, levels = df_barplot$Gene_name %>% rev())

pdf("results/identification_of_candidate_genes/output/top10_shared_DEGs_associated_with_INS_barplot.pdf", width = 3, height = 3)
ggplot(df_barplot, aes(x = log10p, y = Gene_name)) +
  geom_bar(stat = "identity", fill = "orange", width = 0.8) +
  labs(x = expression(paste("-log"[10], "p"))) + 
  theme_bw() +
  theme(axis.title.y = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())
dev.off()


#plot association between DOCK10 and INS expression as in Figure 3D
df_plot_dock10 <- vst_counts[c("ENSG00000254647", res_table[res_table$Gene_name == "DOCK10", "Gene_ID"]), ] %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Sample_ID") %>%
  melt(id.vars = c("Sample_ID", "ENSG00000254647"), variable.name = "Gene_ID") %>%
  merge(annotation, by = "Gene_ID") 

pdf("results/identification_of_candidate_genes/output/DOCK10_INS_association_dotplot.pdf", width = 3.5, height = 3.5)
ggplot(df_plot_dock10) +
  geom_smooth(method = lm, color = "#5658ad", aes(x = value, y = ENSG00000254647)) +
  geom_point(aes(x = value, y = ENSG00000254647)) +
  annotate("text", x = 5.5, y = 20, size = 3,
           label = paste0("Beta = ", res_table[res_table$Gene_name == "DOCK10", "estimate"] %>%
                            round(digits = 2), 
                          ", p = ", res_table[res_table$Gene_name == "DOCK10", "p_value"] %>% 
                            formatC(format = "e", digits = 0))) +
  labs(x = "DOCK10\nvariance stabilized expression", y = "INS\nvariance stabilized expression") +
  theme_bw() +
  theme(panel.grid = element_blank())
dev.off()




#dotplot for top 10 genes associated with INS expression as in Supplementary Figure 3B
df_plot_top10 <- vst_counts[c("ENSG00000254647", res_table[c(1, 3:10), "Gene_ID"]), ] %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Sample_ID") %>%
  melt(id.vars = c("Sample_ID", "ENSG00000254647"), variable.name = "Gene_ID") %>%
  merge(annotation, by = "Gene_ID")

df_plot_top10$Gene_name <- factor(df_plot_top10$Gene_name, levels = res_table[c(1, 3:10), "Gene_name"])

df_min <- df_plot_top10 %>%
  group_by(Gene_name) %>%
  summarize("x" = min(value) + (max(value) - min(value))*0.4)

association_text <- res_table[c(1, 3:10), ] %>%
  merge(df_min, by = "Gene_name")
association_text$Gene_name <- factor(association_text$Gene_name, levels = res_table[c(1, 3:10), "Gene_name"])

pdf("results/identification_of_candidate_genes/output/top10_shared_DEGs_associated_with_INS_dotplot.pdf", width = 8, height = 8)
ggplot(df_plot_top10, aes(x = value, y = ENSG00000254647)) +
  geom_smooth(method = lm, color = "#5658ad") +
  geom_pointassociation_text
  geom_text(data = res_table_text, aes(label = paste0("Beta = ", estimate %>% round(digits = 2), 
                                                      ", p = ", p_value %>% round(digits = 4)), 
                                       x = x, y = 21)) +
  facet_wrap(.~ Gene_name, nrow = 3, scales = "free_x") +
  labs(x = "Variance stabilized expression", y = "INS\nvariance stabilized expression") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
dev.off()



