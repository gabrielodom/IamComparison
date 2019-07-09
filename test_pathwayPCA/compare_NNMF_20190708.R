# Test NNMF Multi-Omics
# Gabriel Odom
# 2019-06-14



######  Description and Overview  #############################################
# We have simulation results from NNMF under the same simulation parameters as
#   the pathwayPCA multi-omics sim. They are
#   - ngroups = 20
#   - p = 1600 (gene expression) & 2400 (DNA methylation)
#   - DE% = 10, 15, 20, 50%
#   - DEsize = 2x, 3x, 4x, 6x
#   - Background DE = 0% (this is a departure from the Pucher paper's 1%)

# I need the genes identified as DE within each run at each design point. Then,
#   I need to group these genes into their original 20 pathways. After that, I
#   perform a Fisher's exact test (example from 50% DE per true positive
#   pathway):
#   
#        In Pathway g     Not in Pathway g     Total
#    DE            n1                   n2       200
#   !DE            n3                   n4      1400
#   Total          80                 1520      1600
#   

# Now that we have a p-value for each pathway, we also need to know the true
#   pathways. For this, we have to call back to the original data, and calculate
#   which pathways are truly DE. Then, finally, we can compare the decisions
#   made with the Fisher's Exact Test with the true pathways. 
# NOTE: the Pucher et al. paper create two such tables, then add each of the
#   counts across data type, yielding
#   
#        In Pathway g     Not in Pathway g     Total
#    DE            N1                   N2       500
#   !DE            N3                   N4      3500
#   Total         200                 3800      4000
#   
#   This equally weights features from both data sets, so a single site is 
#   worth as much as an entire gene -- quite foolish.




######  Setup  ################################################################
library(tidyverse)
library(pathwayPCA)

# The data files all reside in
dataDir <- "results/GeneExp_Met_2DS/"

# Simulated data are in
simDir <- paste0(
  dataDir, "syntheticData/gr_20_perc_10_50_es_20_60_bkgrnd_0/RData/"
)

# Helper functions are
source("test_pathwayPCA/global_functions.R")

# The global parameters are
nRuns_int <- 100L
p.gene.exp <- 1600
p.methylation <- 2400
n.groups <- 20
geneGrpCard_int <- ceiling(p.gene.exp / n.groups)
methylGrpCard_int <- ceiling(p.methylation / n.groups)


# The design component of the file names are
designs_char <- c(
  "20gr_50p_60", "20gr_50p_40", "20gr_50p_30", "20gr_50p_20",
  "20gr_20p_60", "20gr_20p_40", "20gr_20p_30", "20gr_20p_20",
  "20gr_15p_60", "20gr_15p_40", "20gr_15p_30", "20gr_15p_20",
  "20gr_10p_60", "20gr_10p_40", "20gr_10p_30", "20gr_10p_20"
)
desg <- designs_char[1]

results_ls <- lapply(designs_char, function(desg){
  
  
  localParams_num <- 
    str_split(desg, "_") %>% 
    unlist() %>% 
    str_extract("[0-9]*") %>% 
    as.numeric()
  names(localParams_num) <- c("numGroups", "pctDE", "pctEffectSize")
  
  
  
  ######  Pathway Collections  ################################################
  ###  Gene Expression Pathways  ###
  # vector of all gene names
  allGenes_char <- paste0("ge_", seq_len(p.gene.exp))
  
  # list of gene names by pathway
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
  # vector of all gene names
  allSites_char <- paste0("met_", seq_len(p.methylation))
  
  # list of site names by pathway
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
  
  
  ######  Apply over Runs  ####################################################
  
  # NNMF Results are in
  resDir <- paste0(dataDir, "synthetNMF/supervised_k2-ES/files/")
  runs <- paste0("run", seq_len(nRuns_int), "/")
  
  runDir <- runs[1]
  # a <- Sys.time()
  out_ls <- lapply(runs, function(runDir){
    
    ######  Extract New Results  ##############################################
    # We need to know which features NNMF claims are DE.
    
    # NNMF Genes
    geneRes_path <- paste0("NMF_GeneExp_synth_", desg, ".txt")
    NNMFgenes_char <- scan(
      paste0(resDir, runDir, geneRes_path),
      what = character()
    )
    
    # NNMF Methylated Sites
    methylRes_path <- paste0("NMF_Methyl_synth_", desg, ".txt")
    NNMFmethyl_char <- scan(
      paste0(resDir, runDir, methylRes_path),
      what = character()
    )
    
    
    
    ######  Fisher's Exact Test for All Pathways  #############################
    
    # ###  Gene Expression  ###
    # # Genes in pw1 & marked as DE
    # n1 <- length(
    #   intersect(
    #     gene_PC$pathways$path1,
    #     NNMFgenes_char
    #   )
    # )
    # # Genes !(in pw1) & marked as DE
    # n2 <- length(
    #   intersect(
    #     setdiff(allGenes_char, gene_PC$pathways$path1),
    #     NNMFgenes_char
    #   )
    # )
    # # Genes in pw1 & !(marked as DE)
    # n3 <- length(
    #   intersect(
    #     gene_PC$pathways$path1,
    #     setdiff(allGenes_char, NNMFgenes_char)
    #   )
    # )
    # # Genes !(in pw1) & !(marked as DE)
    # n4 <- length(
    #   intersect(
    #     setdiff(allGenes_char, gene_PC$pathways$path1),
    #     setdiff(allGenes_char, NNMFgenes_char)
    #   )
    # )
    # 
    # fisher.test(
    #   matrix(c(n1, n2, n3, n4), ncol = 2, byrow = TRUE),
    #   alternative = "greater"
    # )$p.value
    
    
    ###  Function Time  ###
    # MOVED TO test_pathwayPCA/global_functions.R
    
    # Test
    geCounts_int <- confusion(
      pathway = gene_PC$pathways$path1,
      allFeatures_char = allGenes_char,
      deFeatures_char = NNMFgenes_char
    )
    # pathwayFisherExact(geCounts_int)
    
    metCounts_int <- confusion(
      pathway = methyl_PC$pathways$path1,
      allFeatures_char = allSites_char,
      deFeatures_char = NNMFmethyl_char
    )
    # pathwayFisherExact(metCounts_int)
    # # The sample size is much larger here, so the p-value will be smaller. Is
    # #   it a more "fair" comparison to measure pathway significance based on 
    # #   gene expression alone? Shouldn't a proper multi-omic procedure have
    # #   an integrated pathway p-value?
    # # From their function file drawROC.R, the calculateACC() function simply
    # #   sums the TP, TN, FP, and FN counts for both data sets. Their code is:
    # TP = as.numeric(TP_ge + TP_met)
    # FP = as.numeric(FP_ge + FP_met)
    # TN = as.numeric(TN_ge + TN_met)
    # FN = as.numeric(FN_ge + FN_met)
    # pathwayFisherExact(geCounts_int + metCounts_int)
    
    
    ###  Apply Time  ###
    jointPathways_ls <- mapply(c, 
                               gene_PC$pathways,
                               methyl_PC$pathways,
                               SIMPLIFY = FALSE
    )
    
    NNMFconfusion_ls <- lapply(
      jointPathways_ls,
      confusion,
      allFeatures_char = c(allGenes_char, allSites_char),
      deFeatures_char = c(NNMFgenes_char, NNMFmethyl_char)
    )
    
    NNMF_pVals <- sapply(NNMFconfusion_ls, pathwayFisherExact)
    
    
    
    ######  True DE Pathways  #################################################
    # Now we need to know which features are truly DE
    
    ###  Differentially-expressed genes  ###
    geneDE_path <- paste0("GeneExp_synth_", desg, "_Diff.RData")
    diffGenes_char <- readRDS(paste0(simDir, runDir, geneDE_path))
    
    diffGenes_idx <- 
      diffGenes_char %>% 
      str_replace("ge_", "") %>% 
      as.integer() %>% 
      sort
    
    exprGene_logi <- seq_len(p.gene.exp) %in% diffGenes_idx
    
    # Which gene pathways were expressed?
    signifPathsG_int <- which(
      sapply(
        seq_len(n.groups),
        isGroupDE,
        grpCard_int = geneGrpCard_int,
        DE_logi = exprGene_logi,
        pctDE = localParams_num[2] / 100
      )
    )
    
    
    
    ######  Return  ###########################################################
    data.frame(
      pVals = NNMF_pVals,
      dePath = seq_along(NNMF_pVals) %in% signifPathsG_int
    )
    
    # END for() runs
  })
  # b <- Sys.time()
  # b - a # 1.070069 sec for 15 reps
  
  names(out_ls) <- paste0("run", 1:15)
  out_df <- bind_rows(out_ls, .id = "run")
  # saveRDS(out_df, file = "results/sim_NNMF_20190614/res_100_80.RDS")
  
  out_df
  
  # END for() designs
})

names(results_ls) <- designs_char
results_df <- bind_rows(results_ls, .id = "design")
resultsClean_df <- 
  results_df %>% 
  separate(design, into = c("nGroups", "pctDE", "effSize")) %>% 
  mutate(nGroups = str_replace(nGroups, "gr", "")) %>% 
  mutate(pctDE = str_replace(pctDE, "p", "")) %>% 
  mutate(run = str_replace(run, "run", ""))

write_csv(
  resultsClean_df,
  path = "results/sim_NNMF_20190614/NNMF_all_designs_15runs.csv"
)


# res100_80_df %>% 
#   mutate(claimDE = pVals < 0.05) %>% 
#   mutate(T1err = claimDE & !dePath) %>% 
#   mutate(T2err = (!claimDE) & dePath) %>% 
#   group_by(run) %>% 
#   summarise(
#     pctT1 = mean(T1err),
#     pctT2 = mean(T2err)
#   )
