## Setup ----------------------------------------
### Bioconductor and CRAN libraries used

library(genefilter)
library(geneplotter)
library(topGO)
library(DBI)
library(plyr)


## Function to run GO enrichment analysis--------------

run_GO_ORA <- function(degs, selected_database, annotation, min_node_size = 5){
  
  geneIDs <- rownames(overallBaseMean)
  inUniverse <- geneIDs %in% geneIDs
  inSelection <- geneIDs %in% degs 
  alg <- factor(as.integer(inSelection[inUniverse]))
  names(alg) <- geneIDs[inUniverse]
  
  #GO enrichment analysis
  onts <- c("MF", "BP", "CC")
  tab <- as.list(onts)
  names(tab) <- onts
  for(i in 1:3){
    
    tgd <- new("topGOdata", 
               ontology = onts[i], 
               allGenes = alg, 
               nodeSize = min_node_size,
               annot = annFUN.org, 
               mapping = selected_database, 
               ID = "ENSEMBL")
    
    resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher" )
    resultTopGO.classic <- runTest(tgd, algorithm = "classic", statistic = "Fisher" )
    
    if(length(nodes(graph(tgd))) < 200){
      table_temp <- GenTable( tgd, Fisher.elim = resultTopGO.elim, 
                              Fisher.classic = resultTopGO.classic,
                              orderBy = "Fisher.elim" , topNodes = length(nodes(graph(tgd))))
    }else{
      table_temp <- GenTable( tgd, Fisher.elim = resultTopGO.elim, 
                              Fisher.classic = resultTopGO.classic,
                              orderBy = "Fisher.elim" , topNodes = 200)
    }
    
    #add gene symbols
    my_genes <- table_temp$GO.ID %>%
      genesInTerm(object = tgd)
    
    table_temp$Genes <- my_genes %>%
      map_chr(.f = function(x){intersect(x, sigGenes(tgd)) %>% 
          paste(collapse = ", ") %>% 
          return()})
    
    table_temp$Symbols <- my_genes %>%
      map_chr(.f = function(x){annotation %>%
          subset(.$Gene_ID %in% intersect(x, sigGenes(tgd))) %>%
          .$Gene_name %>% 
          paste(collapse = ", ") %>%
          return()})
    
    table_temp$Fisher.elim <- as.numeric(table_temp$Fisher.elim)
    table_temp$Fisher.classic <- as.numeric(table_temp$Fisher.classic)
    table_temp$Term <- getTermsDefinition(whichTerms = table_temp$GO.ID, onts[i])
    
    tab[[i]] <- table_temp
  }
  
  tab$MF$ont <- "MF"
  tab$BP$ont <- "BP"
  tab$CC$ont <- "CC"
  topGOResults <- rbind.fill(tab)
  
  return(topGOResults)
}




## Function to filter GO results for significance, minimum number of genes, ontology, unique gene sets and maximum number of terms -------------

filter_go_results <- function(go_results, top_num, onto = "ALL"){
  
  go_results <- go_results %>%
    subset(.$Fisher.elim < 0.05) %>%
    subset(.$Significant > 1)
  
  if(onto != "ALL"){
    go_results <- subset(go_results, go_results$ont == onto)
  }
  
  if (nrow(go_results) > 0) {
    go_results <- go_results[order(go_results$Fisher.elim), ]
    exclude_indices <- c()
    unique_gene_sets <- c()
    for(i in 1:length(go_results$GO.ID)){
      if(as.vector(go_results$Genes[i]) %in% unique_gene_sets){
        exclude_indices <- c(exclude_indices, i)
      }else{
        unique_gene_sets <- c(unique_gene_sets, as.vector(go_results$Genes[i]))
      }
    }
    
    if(length(exclude_indices) > 0){
      go_results_filtered <- go_results[-exclude_indices, ]
    }else{
      go_results_filtered <- go_results
    }
    
    if(length(go_results_filtered$GO.ID) > top_num){
      top <- go_results_filtered[1:top_num,]
    }else{
      top <- go_results_filtered
    }
    
    top$Fisher.elim
    top$p <- -1*log10(as.numeric(top$Fisher.elim))
    top$Ratio <- top$Significant / top$Annotated
    
    return(top)
  }
}



## Function to get full GO term names -----------------

getTermsDefinition <- function(whichTerms, ontology, numChar = 1000, multipLines = FALSE) {
  
  qTerms <- paste(paste("'", whichTerms, "'", sep = ""), collapse = ",")
  retVal <- dbGetQuery(GO_dbconn(), paste("SELECT term, go_id FROM go_term WHERE ontology IN",
                                          "('", ontology, "') AND go_id IN (", qTerms, ");", sep = ""))
  
  termsNames <- retVal$term
  names(termsNames) <- retVal$go_id
  
  if(!multipLines)
    shortNames <- paste(substr(termsNames, 1, numChar),
                        ifelse(nchar(termsNames) > numChar, '...', ''), sep = '')
  else
    shortNames <- sapply(termsNames,
                         function(x) {
                           a <- strwrap(x, numChar)
                           return(paste(a, sep = "", collapse = "\\\n"))
                         })
  
  names(shortNames) <- names(termsNames)
  
  #return NAs for the terms that are not found in the DB and make sure the 'names' attribute is as specified
  shortNames <- shortNames[whichTerms]
  names(shortNames) <- whichTerms
  
  return(shortNames)
}




## Function to create a data.frame object that can be used to plot GO terms -------------------------

get_GO_df <- function(GO_tables, direction, score_name, reverse_bool = FALSE, onto = "BP", showTerms = 10, multLines = TRUE, numChar = 60, return_level_list = FALSE){
  
  #get top terms 
  table <- GO_tables
  score_index <- which(names(table[[1]]) == score_name)
  i <- 1
  df_terms <- c()
  df_terms_list <- c()
  for(i in 1:length(GO_tables)) {
    colnames(table[[i]])[score_index] <- "Scores"
    
    if(onto != "ALL"){
      table[[i]] <- table[[i]] %>%
        subset(.$ont == onto)
    }
    
    #Order data frame by P-value
    idx <- order(table[[i]]$Scores, decreasing = FALSE)
    table[[i]] <- table[[i]][idx, ]
    
    ifelse(is.numeric(showTerms), 
           table[[i]] <- table[[i]][1:showTerms, ],
           table[[i]] <- table[[i]] %>% subset(.$GO.ID %in% showTerms))
    df_terms <- union(df_terms, table[[i]]$GO.ID)
    df_terms_list[[i]] <- table[[i]]$GO.ID
  }
  names(df_terms_list) <- direction
  
  df <- data.frame()
  for(i in 1:length(GO_tables)) {
    df_temp <- GO_tables[[i]][GO_tables[[i]]$GO.ID %in% df_terms, ]
    df_temp <- df_temp %>%
      add_column(Direction = i)
    
    df <- rbind(df, df_temp)
  }
  
  if(reverse_bool){
    my_levels <- table %>%
      lapply(function(x){x$Term %>% rev()}) %>%
      c(recursive = TRUE) %>%
      unique()
  }else if(!reverse_bool){
    my_levels <- table %>%
      lapply(function(x){x$Term}) %>%
      c(recursive = TRUE) %>%
      unique() %>%
      rev()
  }
  
  df$Term <- factor(df$Term, levels = my_levels)
  colnames(df)[7] <- "Scores"
  df$Direction <- factor(df$Direction)
  levels(df$Direction) <- direction
  
  #reduce length of term names 
  termsNames <- levels(df$Term)
  
  if(multLines == FALSE) {
    shortNames <- paste(substr(as.character(termsNames), 1, numChar),
                        ifelse(nchar(as.character(termsNames)) > numChar, '...', ''), sep = '')
  } else {
    shortNames <- sapply(as.character(termsNames),
                         function(x) {
                           a <- strwrap(x, numChar)
                           return(paste(a, sep = "", collapse = "\n"))
                         })
  }
  
  levels(df$Term) <- shortNames
  
  df <- df %>% 
    add_column(changed_scores = -1*log10(df$Scores))
  
  if(return_level_list){
    return(list(df, df_terms_list))
  }else if(!return_level_list){
    return(df) 
  }
}



## Function to create a data.frame object that can be used to plot RGD terms -------------------------

get_KEGG_df <- function(GO_tables, direction, score_name, reverse_bool = FALSE, showTerms = 10, multLines = TRUE, numChar = 60, return_level_list = FALSE){
  
  #get top terms 
  table <- GO_tables
  score_index <- which(names(table[[1]]) == score_name)
  i <- 1
  df_terms_list <- 
  df_terms <- c()
  for(i in 1:length(GO_tables)) {
    colnames(table[[i]])[score_index] <- "Scores"
    
    # Order data frame by P-value
    idx <- order(table[[i]]$Scores, decreasing = FALSE)
    table[[i]] <- table[[i]][idx, ]
    
    ifelse(is.numeric(showTerms), 
           table[[i]] <- table[[i]][1:showTerms, ],
           table[[i]] <- table[[i]] %>% subset(.$ID %in% showTerms))
    df_terms <- union(df_terms, table[[i]]$ID)
    df_terms_list[[i]] <- table[[i]]$ID
  }
  names(df_terms_list) <- direction
  
  #get all data
  df <- data.frame()
  for(i in 1:length(GO_tables)) {
    df_temp <- GO_tables[[i]][GO_tables[[i]]$ID %in% df_terms, ]
    
    #add group
    df_temp <- df_temp %>%
      add_column(Direction = i)
    
    df <- rbind(df, df_temp)
  }
  
  if(reverse_bool){
    my_levels <- table %>%
      lapply(function(x){x$Description %>% rev()}) %>%
      c(recursive = TRUE) %>%
      unique()
  }else if(!reverse_bool){
    my_levels <- table %>%
      lapply(function(x){x$Description}) %>%
      c(recursive = TRUE) %>%
      unique() %>%
      rev()
  }
  
  df$Description <- factor(df$Description, levels = my_levels)
  colnames(df)[7] <- "Scores"
  df$Direction <- factor(df$Direction)
  levels(df$Direction) <- direction
  
  #reduce length of term names 
  termsNames <- levels(df$Description)
  
  if(multLines == FALSE) {
    shortNames <- paste(substr(as.character(termsNames), 1, numChar),
                        ifelse(nchar(as.character(termsNames)) > numChar, '...', ''), sep = '')
  } else {
    shortNames <- sapply(as.character(termsNames),
                         function(x) {
                           a <- strwrap(x, numChar)
                           return(paste(a, sep = "", collapse = "\n"))
                         })
  }
  
  levels(df$Description) <- shortNames
  
  #add -log10 pvalue and set sizes of terms
  df$setSize <- df$BgRatio %>%
    str_split("/") %>%
    do.call(rbind, .) %>%
    as.data.frame() %>%
    .[, 1] %>%
    as.numeric()
  df <- df %>% 
    add_column("changed_scores" = -1*log10(df$Scores),
               "Ratio" = .$Count/.$setSize)
  
  if(return_level_list){
    return(list(df, df_terms_list))
  }else if(!return_level_list){
    return(df) 
  }
}

