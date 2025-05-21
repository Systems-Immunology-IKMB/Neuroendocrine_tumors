## Setup ----------------------------------------
### Bioconductor and CRAN libraries used

library(tidyverse)
library(tximport)
library(DESeq2)

setwd("C:/Projects/Neuroendocrine_tumors")



## Principal component analysis for surgical samples ------------------------------------------------------------

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


#run PCA
pca <- DESeqDataSetFromTximport(txi, colData = metadata, design = ~ 1) %>%
  .[rowSums(counts(.) == 0) < 0.8*nrow(metadata), ] %>%
  estimateSizeFactors() %>%
  .[, metadata$Sample_ID] %>%
  vst() %>%
  assay() %>%
  t() %>%
  prcomp()

pc_variance <- pca %>%
  summary() %>%
  .$importance %>%
  .["Proportion of Variance", ] %>% 
  `*`(100) %>%
  round(2)

df_pca <- cbind(metadata, pca$x)


#plot results for PC1 and PC2 as in Figure 1E
pdf("results/PCA/output/PCA_1_2_colorPart_shapeGroup_surgical_samples.pdf", width = 7, height = 5)
ggplot(df_pca, aes(x = PC1, y = PC2)) +
  geom_point(size = 1.9, aes(color = Part, shape = Group)) +
  scale_shape_manual(values = c("Insulinoma" = 3, "NF-NET" = 16)) +
  scale_color_manual(values = c("Normal" = "#66C2A5", "Tumor" = "#FFAA42", "Liver metastasis" = "#DC6789")) +
  labs(x = paste0("PC1 (", pc_variance[1] %>% round(2), "% variance)"),
       y = paste0("PC2 (", pc_variance[2]%>% round(2), "% variance)"),
       color = "Surgical sample",
       shape = "Group") +
  theme_bw() +
  theme(panel.border = element_rect(color = "black"),
        panel.grid = element_blank())
dev.off()







## Principal component analysis for organoid samples ------------------------------------------------------------

#load metadata
metadata <- read.csv("info_files/Sequenced_sample_info_organoid.csv")
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

files <- list.files(path = "count_files_organoid_samples/salmon_files", full.names = T, pattern = "quant.sf", recursive = TRUE)
names(files) <- str_replace(files, "count_files_organoid_samples/salmon_files/", "") %>% 
  str_replace("_quant.sf", "") 
files <- files[metadata$Sample_ID]
txi <- tximport(files, type = "salmon", tx2gene = tx2gene[, 1:2], countsFromAbundance = "no")

all(colnames(txi$counts) == metadata$Sample_ID) #TRUE


#run PCA
pca <- DESeqDataSetFromTximport(txi, colData = metadata, design = ~ 1) %>%
  .[rowSums(counts(.) == 0) < 0.8*nrow(metadata), ] %>%
  estimateSizeFactors() %>%
  .[, metadata$Sample_ID] %>%
  vst() %>%
  assay() %>%
  t() %>%
  prcomp()

pc_variance <- pca %>%
  summary() %>%
  .$importance %>%
  .["Proportion of Variance", ] %>% 
  `*`(100) %>%
  round(2)

df_pca <- cbind(metadata, pca$x)


#plot results for PC1 and PC2 as in Figure 1F
pdf("results/PCA/output/PCA_1_2_colorPart_shapeGroup_organoid_samples.pdf", width = 7, height = 5)
ggplot(df_pca, aes(x = PC1, y = PC2)) +
  geom_point(size = 1.9, aes(color = Part, shape = Group)) +
  scale_shape_manual(values = c("Insulinoma" = 3, "NF-NET" = 16)) +
  scale_color_manual(values = c("Normal" = "#66C2A5", "Tumor" = "#FFAA42", "Liver metastasis" = "#DC6789")) +
  labs(x = paste0("PC1 (", pc_variance[1] %>% round(2), "% variance)"),
       y = paste0("PC2 (", pc_variance[2]%>% round(2), "% variance)"),
       color = "Surgical sample",
       shape = "Group") +
  theme_bw() +
  theme(panel.border = element_rect(color = "black"),
        panel.grid = element_blank())
dev.off()

