## Setup ----------------------------------------
### Bioconductor and CRAN libraries used

library(tidyverse)
library(reshape2)
library(tximport)
library(DESeq2)
library(genefilter)
library(geneplotter)
library(topGO)
library(DBI)
library(plyr)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)

setwd("C:/Projects/Neuroendocrine_tumors")
source("C:/Projects/Neuroendocrine_tumors/scripts/enrichment_analysis_helper_functions.R")


##Loading data ------------------------------------------------------------

tx2gene <-  read.table("info_files/salmon_tx2gene.tsv")
colnames(tx2gene) <- c("Transcript_ID", "Gene_ID", "Gene_name")
annotation <- tx2gene[, 2:3] %>% unique()
rownames(annotation) <- annotation$Gene_ID





## Functional enrichment analysis using RGD --------------------------------------------

#import RGD data (Human Gene-to-Pathway Annotations File downloaded from https://rgd.mcw.edu/wg/home/pathway2/)
table_rgd_pathways <- read.table("info_files/RGD_gene2pathway.txt", header = TRUE, sep = '\t', fill = TRUE, quote = "") %>%
  merge(annotation, by.x = "OBJECT_SYMBOL", by.y = "Gene_name")

df_rgd <- table_rgd_pathways[, c("TERM_ACC_ID", "Gene_ID")]
colnames(df_rgd) <- c("term", "gene")

rgd_term_annotation <- table_rgd_pathways[, c("TERM_ACC_ID", "TERM_NAME")] %>%
  unique()


#get differentially expressed genes
folder <- paste0("results/DGE_analysis/output/")
res_table <- read.table(paste0(folder, "DEres_surgical_samples_Insulinoma_Normal_vs_Tumor.txt"), header = TRUE, sep = "\t")
res_table <- res_table %>%
  rownames_to_column("Gene_ID")
res_table$Gene_ID <- res_table$Gene_ID %>%
  substr(1, 18)
rownames(res_table) <- res_table$Gene_ID

degs_tab <- res_table %>%
  subset(.$padj < 0.05) %>%
  subset(.$baseMean > 10) %>%
  subset(abs(.$log2FoldChange) > 1)


#run enrichment analysis for up- and down-regulated DEGs
directions <- c("up", "down")
for(direction in  directions){
  
  if(direction == "up"){
    degs <- degs_tab %>%
      subset(.$log2FoldChange > 0) %>%
      rownames()
  }else if(direction == "down"){
    degs <- degs_tab %>%
      subset(.$log2FoldChange < 0) %>%
      rownames()
  }

  geneIDs <- rownames(res_table)
  inSelection <-  geneIDs %>% 
    subset(geneIDs %in% degs)
  
  enrich_data <- enricher(gene = inSelection,
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          universe = geneIDs,
                          qvalueCutoff = 0.05,
                          TERM2GENE = df_rgd)
  
  res <- enrich_data@result
  
  #add term info
  res <- res %>%
    merge(rgd_term_annotation, by.x = "ID", by.y = "TERM_ACC_ID", sort = FALSE) 
  
  my_genes <- res$geneID %>%
    sapply(FUN = function(x){x %>% str_split(pattern = "/") %>% return()})
  names(my_genes) <- res$ID
  
  res$geneID <- my_genes %>%
    map_chr(.f = function(x){x %>% 
        paste(collapse = ", ") %>% 
        return()})
  
  res$Symbols <- my_genes %>%
    map_chr(.f = function(my_genes){annotation %>%
        subset(.$Gene_ID %in% my_genes) %>%
        .$Gene_name %>% 
        paste(collapse = ", ") %>%
        return()})
  
  #filter for significant results and export
  res_sig <- res %>%
    subset(.$qvalue < 0.05)
  
  file_name <- paste0(folder, "/ORA_RGD_Insulinoma_", direction, "_degs.csv")
  write.csv(res_sig, file_name, row.names = FALSE)
}






## Functional enrichment analysis using GO terms --------------------------------------------

#get differentially expressed genes
folder <- paste0("results/DGE_analysis/output/")
res_table <- read.table(paste0(folder, "DEres_surgical_samples_Insulinoma_Normal_vs_Tumor.txt"), header = TRUE, sep = "\t")
res_table <- res_table %>%
  rownames_to_column("Gene_ID")
res_table$Gene_ID <- res_table$Gene_ID %>%
  substr(1, 15)

degs_tab <- res_table %>%
  subset(.$padj < 0.05) %>%
  subset(.$baseMean > 10) %>%
  subset(abs(.$log2FoldChange) > 1)


#run enrichment analysis for up- and down-regulated DEGs
directions <- c("up", "down")
for(direction in  directions){
  
  if(direction == "up"){
    degs <- degs_tab %>%
      subset(.$log2FoldChange > 0) %>%
      .$Gene_ID %>% 
      unique()
  }else if(direction == "down"){
    degs <- degs_tab %>%
      subset(.$log2FoldChange < 0) %>%
      .$Gene_ID %>% 
      unique()
  }
  
  overallBaseMean <- as.matrix(res_table[, "baseMean", drop = F])
  colnames(overallBaseMean) <- "mean_expression"
  rownames(overallBaseMean) <- res_table$Gene_ID %>% 
    str_sub(1, 15)
  
  topGOResults <- run_GO_ORA(degs, 
                             overallBaseMean, 
                             selected_database = "org.Hs.eg.db", 
                             back_num = 10,
                             background = FALSE,
                             annotation = annotation)
  
  go_results_filtered <- filter_go_results(topGOResults, top_num = 600, onto = "BP")
  
  file_filtered <- paste0(folder, "/ORA_GO_Insulinoma_", direction, "_degs.csv")
  write.csv(go_results_filtered, file_filtered, row.names = FALSE)
}





## ORA dotplot as in Figure 2C -------------------------------------------------

directions <- c("up", "down")

#get GO terms
GO_tables <- list()
section <- c()
for(direction in directions){
  
  folder <- paste0("results/DGE_analysis/output/")
  file_filtered <- paste0(folder, "/ORA_GO_Insulinoma_", direction, "_degs.csv")
  my_table <- read_csv(file_filtered) 
  my_table <- my_table %>% 
    subset(.$Fisher.elim < 0.05) 
  
  GO_tables <- append(GO_tables,
                      list(my_table))
  if(nrow(my_table > 0)){section <- c(section, paste0("Insulinoma_vs_normal,\n", direction, "regulated"))}
}

df_plot_GO <- get_GO_df(GO_tables = GO_tables,
                        direction = section,
                        reverse_bool = FALSE,
                        score_name = "Fisher.elim",
                        return_level_list = TRUE,
                        showTerms = 5,
                        multLines = TRUE,
                        numChar = 60)

#get RGD terms
RGD_tables <- list()
section <- c()
for(direction in directions){
  
  folder <- paste0("results/DGE_analysis/output/")
  file_filtered <- paste0(folder, "/ORA_RGD_Insulinoma_", direction, "_degs.csv")
  my_table <- read_csv(file_filtered) 
  my_table <- my_table %>% 
    subset(.$qvalue < 0.05) 
  my_table$Description <- my_table$TERM_NAME
  
  RGD_tables <- append(RGD_tables,
                      list(my_table))
  if(nrow(my_table > 0)){section <- c(section, paste0(groups_names[i], "_vs_normal,\n", direction, "regulated"))}
}

df_plot_RGD <- get_KEGG_df(GO_tables = RGD_tables,
                           direction = section,
                           reverse_bool = FALSE,
                           score_name = "qvalue",
                           return_level_list = TRUE,
                           showTerms = 5,
                           multLines = TRUE,
                           numChar = 80)

#merge tables
df_plot_GO_terms <- df_plot_GO[[1]] %>%
  dplyr::rename("ID" = "GO.ID",
                "Term" = "Term",
                "Annotated" = "Annotated",
                "Significant" = "Significant",
                "Ratio" = "Ratio",
                "Scores" = "Scores",
                "Genes" = "Genes",
                "Symbols" = "Symbols",
                "Direction" = "Direction",
                "log10scores" = "changed_scores") %>%
  .[, c("ID", "Term", "Annotated", "Significant", "Ratio", "Scores", "Genes", "Symbols", "Direction", "log10scores")]
df_plot_GO_terms$Term <- paste("GO:", df_plot_GO_terms$Term)

df_plot_RGD_terms <- df_plot_RGD[[1]] %>%
  dplyr::rename("ID" = "ID",
                "Term" = "TERM_NAME",
                "Annotated" = "setSize",
                "Significant" = "Count",
                "Ratio" = "Ratio",
                "Scores" = "Scores",
                "Genes" = "geneID",
                "Symbols" = "Symbols",
                "Direction" = "Direction",
                "log10scores" = "changed_scores") %>%
  .[, c("ID", "Term", "Annotated", "Significant", "Ratio", "Scores", "Genes", "Symbols", "Direction", "log10scores")]
df_plot_RGD_terms$Term <- paste("RGD:", df_plot_RGD_terms$Term)

df_plot_merged <- df_plot_GO_terms %>%
  rbind(df_plot_RGD_terms)

#add term order
section_temp <- section[1]
my_levels <- c()
for(section_temp in section){
  
  my_levels_temp <- c(df_plot_GO[[2]][[section_temp]], df_plot_RGD[[2]][[section_temp]])
  df_terms_temp <- df_plot_merged %>%
    subset(.$Direction == section_temp) %>%
    subset(.$ID %in% my_levels_temp) %>%
    .[order(.$Scores, decreasing = FALSE), ] 
  my_levels <- c(my_levels, df_terms_temp$Term)
}

df_plot_merged$Term <- factor(df_plot_merged$Term, levels = my_levels %>% unique() %>% rev())

#select terms
df_plot_merged <- df_plot_merged %>%
  subset(!.$ID %in% c("PW:0000101", "PW:0001066"))


pdf(paste0("results/DGE_analysis/output/dotplot_ORA_Insulinoma_DEGs.pdf"), height = 8, width = 8)
ggplot(df_plot_merged, mapping = aes(x = Direction, y = Term, size = Ratio, color = log10scores)) +
  geom_point() +
  labs(x = "Normal vs Tumor", y = "GO terms") +
  scale_colour_gradient(high = "#990000", low = "#FF9999") +
  scale_x_discrete(labels = c("Insulinoma_vs_normal,\nupregulated" = "Upregulated in Insulinoma", 
                              "Insulinoma_vs_normal,\ndownregulated" = "Downregulated in Insulinoma")) +
  theme_bw() + 
  theme(axis.text.y = element_text(hjust = 1, size=13, color = "black"), 
        axis.text.x = element_text(size=13, color = "black", angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13)) +
  labs(color = expression(paste("-log"[10], "q")), 
       size = paste0("Gene ratio"))
dev.off()
