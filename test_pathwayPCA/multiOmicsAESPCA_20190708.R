# Simulation Run: pathwayPCA
# Gabriel Odom
# 2019-07-08

# This simulation run is a slimmed-down version of multiOmicsRes_20190603.R.
#   Within each design point and run, this simulation script will create the
#   Omics* object for each data type, calculate the pathways (and features)
#   that are truly differentially-expressed, and calculate the AESPCA results
#   within each design point. We will save the Omics* object, the true pathways,
#   and the AESPCA results for both gene expression and methylation data sets.
# We analyse the results in multiOmicsAESPCAresults_20190708.R.


library(tidyverse)
library(pathwayPCA)



######  Global Functions  #####################################################
source("test_pathwayPCA/global_functions.R")


######  Global Parameters  ####################################################
pheno_char <- rep(c("tumour", "normal"), each = 100)
p.gene.exp <- 1600
p.methylation <- 2400
n.groups <- 20
geneGrpCard_int <- ceiling(p.gene.exp / n.groups)
methylGrpCard_int <- ceiling(p.methylation / n.groups)

whereDir <- 
  "results/GeneExp_Met_2DS/syntheticData/gr_20_perc_10_50_es_20_60_bkgrnd_0/RData/"
runDirs_char <- paste0(list.files(whereDir), "/")


###  Gene Expression Pathways  ###
genePaths_ls <- lapply(
  seq_len(n.groups),
  makePath,
  grpCard_int = geneGrpCard_int,
  pathType = "ge_"
)
names(genePaths_ls) <- paste0("path", seq_along(genePaths_ls))

gene_PC <- CreatePathwayCollection(
  sets_ls = genePaths_ls,
  TERMS = names(genePaths_ls)
)

rm(genePaths_ls)


###  Methlyation Pathways  ###
methylPaths_ls <- lapply(
  seq_len(n.groups),
  makePath,
  grpCard_int = methylGrpCard_int,
  pathType = "met_"
)
names(methylPaths_ls) <- paste0("path", seq_along(methylPaths_ls))

methyl_PC <- CreatePathwayCollection(
  sets_ls = methylPaths_ls,
  TERMS = names(methylPaths_ls)
)

rm(methylPaths_ls)


######  Design Parameters and Loop  ###########################################
designs_char <- c(
  "20gr_50p_60", "20gr_50p_40", "20gr_50p_30", "20gr_50p_20",
  "20gr_20p_60", "20gr_20p_40", "20gr_20p_30", "20gr_20p_20",
  "20gr_15p_60", "20gr_15p_40", "20gr_15p_30", "20gr_15p_20",
  "20gr_10p_60", "20gr_10p_40", "20gr_10p_30", "20gr_10p_20"
)

a <- Sys.time()

# desg <- 1
# designs_char[desg]
# # 20 groups; 50 percent of features DE in perturbed pathways; 6x effect size
for(desg in seq_along(designs_char)[-1]){
  
  localParams_num <- 
    str_split(designs_char[desg], "_") %>% 
    unlist() %>% 
    str_extract("[0-9]*") %>% 
    as.numeric()
  names(localParams_num) <- c("numGroups", "pctDE", "pctEffectSize")
  
  
  ###  Filenames  ###
  geneDE_path <- paste0("GeneExp_synth_", designs_char[desg], "_Diff.RData")
  geneData_path <- paste0("GeneExp_synth_", designs_char[desg], ".RData")
  diffMethyl_path <- paste0("Methyl_synth_", designs_char[desg], "_Diff.RData")
  methylData_path <- paste0("Methyl_synth_", designs_char[desg], ".RData")
  
  
  
  ######  Runs Loop  ##########################################################
  
  # runDir <- runDirs_char[1]
  
  for(runDir in runDirs_char){
    
    ######  Differentially-Expressed Genes  ###################################
    diffGenes_char <- readRDS(
      paste0(whereDir, runDir, geneDE_path)
    )
    
    diffGenes_idx <- 
      diffGenes_char %>% 
      str_replace("ge_", "") %>% 
      as.integer() %>% 
      sort
    
    exprGene_logi <- seq_len(p.gene.exp) %in% diffGenes_idx
    
    # Which pathways were expressed?
    signifPaths_int <- which(
      sapply(
        seq_len(n.groups),
        isGroupDE,
        grpCard_int = geneGrpCard_int,
        DE_logi = exprGene_logi,
        pctDE = localParams_num[2] / 100
      )
    )
    
    rm(diffGenes_idx, exprGene_logi)
    
    
    
    ######  Genes  ############################################################
    gene_mat <- readRDS(
      paste0(whereDir, runDir, geneData_path)
    )
    
    tGene_df <- data.frame(
      SampleIDs = colnames(gene_mat),
      Type = pheno_char,
      t(gene_mat),
      row.names = NULL,
      stringsAsFactors = FALSE
    )
    
    
    ######  Gene OmicsCateg  ##################################################
    gene_Omics <- CreateOmics(
      assayData_df = tGene_df[, -2],
      pathwayCollection_ls = gene_PC,
      response = tGene_df[, 1:2],
      respType = "categ"
    )
    
    gene_aespcOut <- AESPCA_pVals(
      gene_Omics,
      parallel = TRUE,
      numCores = 5,
      adjustment = "BH"
    )
    
    # Time in seconds
    # 29.47 for 2 cores; 22.23 for 4; 21.62 for 5; 24.28 for 6; 26.48 for 8
    
    rm(gene_mat, tGene_df)
    
    
    
    ######  Differentially-Methylated Genes  ##################################
    
    diffMethyl_char <- readRDS(
      paste0(whereDir, runDir, diffMethyl_path)
    )
    
    
    
    ######  Methylation  ######################################################
    methyl_mat <- readRDS(
      paste0(whereDir, runDir, methylData_path)
    )
    
    tMethyl_df <- data.frame(
      SampleIDs = colnames(methyl_mat),
      Type = pheno_char,
      t(methyl_mat),
      row.names = NULL,
      stringsAsFactors = FALSE
    )
    
    
    ######  Methylation OmicsCateg  ###########################################
    methyl_Omics <- CreateOmics(
      assayData_df = tMethyl_df[, -2],
      pathwayCollection_ls = methyl_PC,
      response = tMethyl_df[, 1:2],
      respType = "categ"
    )
    
    methyl_aespcOut <- AESPCA_pVals(
      methyl_Omics,
      parallel = TRUE,
      numCores = 5,
      adjustment = "BH"
    )
    # 23.00 seconds over 5 cores
    
    rm(methyl_mat, tMethyl_df)
    
    
    
    ######  Simulation Experiment Output  #####################################
    simResults_ls <- list(
      fileNames = c(geneData_path, methylData_path),
      parameters = localParams_num,
      truePaths = signifPaths_int,
      geneOmics = gene_Omics,
      geneResults = gene_aespcOut,
      trueGenes = diffGenes_char,
      methylOmics = methyl_Omics,
      methylResults = methyl_aespcOut,
      trueMethyl = diffMethyl_char
    )
    
    out_path <- "results/sim_20190708/"
    outFile_char <- paste(
      "multiOmics",
      paste(localParams_num, collapse = "_"),
      str_extract(runDir, "(\\d)+"), 
      sep = "_"
    )
    
    saveRDS(simResults_ls, file = paste0(out_path, outFile_char, ".RDS"))
    
    
    
    closeAllConnections()
    
    # END for() Runs
    
  }
  
  # END for designs
  
}

Sys.time() - a

# For 1 design point iterating over 100 simulation runs: 1.254376 hrs
# For 16 design point iterating over 100 simulation runs: 18.95096 hours
