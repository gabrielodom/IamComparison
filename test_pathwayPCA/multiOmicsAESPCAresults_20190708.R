
library(tidyverse)
library(pathwayPCA)

######  Global Functions & Parameters  ########################################
source("test_pathwayPCA/global_functions.R")

in_path <- "results/sim_20190708/"
aespcaFiles_char <- list.files(in_path)

n.groups <- 20
pValues_num <- c(0.001, seq(0.005, 0.395, by = 0.005), 0.499)

# file <- aespcaFiles_char[1]
a <- Sys.time()
for (file in aespcaFiles_char) {
  
  
  ######  Data Setup  #########################################################
  aespcaResults_ls <- readRDS(file = paste0(in_path, file))
  
  
  signifPaths_int <- aespcaResults_ls$truePaths
  gene_Omics <- aespcaResults_ls$geneOmics
  gene_aespcOut <- aespcaResults_ls$geneResults
  methyl_Omics <- aespcaResults_ls$methylOmics
  methyl_aespcOut <- aespcaResults_ls$methylResults
  
  ######  Sensitivity and Specificity  ######################################
  geneSnS_df <- lapply(pValues_num, function(pVal){
    
    decisionTable(
      PathpVals_df = getPathpVals(gene_aespcOut),
      truePaths = names(signifPaths_int),
      type = "GeneExp",
      pValThresh = pVal
    )
    
  }) %>% 
    bind_rows
  
  methylSnS_df <- lapply(pValues_num, function(pVal){
    
    decisionTable(
      PathpVals_df = getPathpVals(methyl_aespcOut),
      truePaths = names(signifPaths_int),
      type = "Methylation",
      pValThresh = pVal
    )
    
  }) %>% 
    bind_rows
  
  multiOmicsSnS <-
    rbind(geneSnS_df, methylSnS_df) %>% 
    arrange(pValue) %>% 
    mutate(TPR = TP / (TP + FN)) %>% 
    mutate(FPR = FP / (FP + TN)) %>% 
    mutate(ACC = (TP + TN) / n.groups) %>% 
    mutate(F1 = 2 * TP / (2 * TP + FP + FN)) %>% 
    mutate(
      MCC = (TP * TN - FP * FN) / 
        sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    )
  
  rm(geneSnS_df, methylSnS_df)
  
  
  
  ######  Global p-Value  ###################################################
  
  pathSignif(
    getPathwayCollection(gene_Omics)$TERMS[1],
    getResponse(gene_Omics),
    omicsOut1 = gene_aespcOut,
    omicsOut2 = methyl_aespcOut
  )
  
  multiOmicspVals_df <- lapply(
    getPathwayCollection(gene_Omics)$TERMS,
    pathSignif,
    resp = getResponse(gene_Omics),
    omicsOut1 = gene_aespcOut,
    omicsOut2 = methyl_aespcOut
  ) %>% 
    bind_rows(.id = "pathway") %>% 
    mutate(DEpath = pathway %in% names(signifPaths_int)) %>% 
    arrange(global)
  
  
  
  ######  Simulation Experiment Output  #########################################
  simResults_ls <- list(
    fileName = file,
    parameters = aespcaResults_ls$parameters,
    classification = multiOmicsSnS,
    multiOmicsSignif = multiOmicspVals_df
  )
  
  outFile_char <- paste0(
    "multiOmicsPerformance",
    gsub("multiOmics", "", file)
  )
  
  saveRDS(simResults_ls, file = paste0(in_path, outFile_char))
  
  
  closeAllConnections()
  # 2.895405 sec for 81 p-value levels (from 0-0.4)
  
}
Sys.time() - a