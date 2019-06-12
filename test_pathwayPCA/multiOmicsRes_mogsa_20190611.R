# Simulation Run: mogsa
# Gabriel Odom
# 2019-06-12

# We have already applied pathwayPCA to the synthetic data, analysed pathwayPCA
#   performance, applied sCCA to the synthetic data, and analysed sCCA 
#   performance. Now we need to apply mogsa to the synthetic data.
# mogsa vignette:
# https://bioconductor.org/packages/release/bioc/vignettes/mogsa/inst/doc/mogsa-knitr.pdf

# We have simulation results from sCCA under the same simulation parameters as
#   the pathwayPCA multi-omics sim. They are
#   - ngroups = 20
#   - p = 1600 (gene expression) & 2400 (DNA methylation)
#   - n = 100 (case) & 100 (control)
#   - DE% = 10, 20, 50, 100%
#   - DEsize = 2x, 4x, 8x
#   - Background DE = 0% (this is a departure from the Pucher paper's 1%)

# mogsa finds subject-specific p-values within each pathway. Therefore, if we
#   have 20 pathways and 200 subjects, then mogsa() will return a 20 x 200
#   matrix of p-values. In their vignette, they find the 20 "most significant"
#   pathways by
#   1. creating a 20 x 200 indicator matrix of significance (at 0.01, but we 
#      will use 0.05), 
#   2. summing the rows of this matrix to yield the number of subjects with
#      significant enrichment in each pathway, and
#   3. sorting to the pathways with the 20 largest counts.
# This is not a pathway p-value. However, I believe that, because the subjects
#   are independent, the count of significant subjects within each pathway
#   should follow a Binomial(200, 0.05) distribution under H0. We can find the
#   probability of observing the number of significant subjects or more via the
#   CDF of this Binomial distribution. Steven confirmed that this approach is
#   appropriate.

# BiocManager::install("mogsa")

library(tidyverse)
library(mogsa)

# mogsa data format:
# data(NCI60_4arrays)
# list of p x n data frames; symbols are row names; sample IDs are column names



######  Global Functions  #####################################################
source("test_pathwayPCA/global_functions.R")


######  Global Parameters  ####################################################
n.samples <- 200
p.gene.exp <- 1600
p.methylation <- 2400
n.groups <- 20
geneGrpCard_int <- ceiling(p.gene.exp / n.groups)
methylGrpCard_int <- ceiling(p.methylation / n.groups)
signif.level <- 0.05

whereDir <- 
  "results/GeneExp_Met_2DS/syntheticData/gr_20_perc_10_100_es_20_80_bkgrnd_0/RData/"
runDirs_char <- paste0(list.files(whereDir), "/")


###  Shared Pathway Collection  ###
# Gene Expression Pathways
genePaths_ls <- lapply(
  seq_len(n.groups),
  makePath,
  grpCard_int = geneGrpCard_int,
  pathType = "ge_"
)

# Methlyation Pathways
methylPaths_ls <- lapply(
  seq_len(n.groups),
  makePath,
  grpCard_int = methylGrpCard_int,
  pathType = "met_"
)

# Combine
mogsaPathways_ls <- mapply(c, genePaths_ls, methylPaths_ls, SIMPLIFY = FALSE)
names(mogsaPathways_ls) <- paste0("path", seq_along(mogsaPathways_ls))

rm(genePaths_ls, methylPaths_ls)



######  Design Parameters and Loop  ###########################################
designs_char <- c(
  "20gr_100p_80", "20gr_100p_40", "20gr_100p_20",
  "20gr_50p_80", "20gr_50p_40", "20gr_50p_20",
  "20gr_20p_80", "20gr_20p_40", "20gr_20p_20",
  "20gr_10p_80", "20gr_10p_40", "20gr_10p_20"
)
# desg <- 1
# designs_char[desg]
# # 20 groups; 100 percent of features DE in perturbed pathways; 80% effect size

a <- Sys.time()

for(desg in seq_along(designs_char)){
  
  localParams_num <- 
    str_split(designs_char[desg], "_") %>% 
    unlist() %>% 
    str_extract("[0-9]*") %>% 
    as.numeric()
  names(localParams_num) <- c("numGroups", "pctDE", "pctEffectSize")
  
  
  ###  Filenames  ###
  geneData_path <- paste0("GeneExp_synth_", designs_char[desg], ".RData")
  methylData_path <- paste0("Methyl_synth_", designs_char[desg], ".RData")
  
  geneDE_path <- paste0("GeneExp_synth_", designs_char[desg], "_Diff.RData")
  methylDE_path <- paste0("Methyl_synth_", designs_char[desg], "_Diff.RData")
  
  
  
  ######  Runs Loop  ##########################################################
  
  # runDir <- runDirs_char[1]
  
  # for(runDir in runDirs_char){
  a <- Sys.time()
  simResults_ls <- lapply(runDirs_char, function(runDir){
    
    ######  Data Array List  ##################################################
    # Genes
    gene_mat <- readRDS(
      paste0(whereDir, runDir, geneData_path)
    )
    
    gene_df <- data.frame(
      gene_mat
    )
    
    
    # Methylation Sites
    methyl_mat <- readRDS(
      paste0(whereDir, runDir, methylData_path)
    )
    
    methyl_df <- data.frame(
      methyl_mat
    )
    
    
    # Combined List
    multiOmics_ls <- list(
      geneExpr  = gene_df,
      dnaMethyl = methyl_df
    )
    
    rm(gene_mat, gene_df, methyl_mat, methyl_df)
    
    
    
    ######  mogsa  ############################################################
    # Pathway Collection to Gene Set Annotation Matrix
    geneSetMat_ls <- prepSupMoa(
      multiOmics_ls,
      geneSets = mogsaPathways_ls,
      minMatch = 1
    )
    
    # mogsa run
    # According to the help files, nf, proc.row, and w.data must be supplied.
    #   The nf parameter "should be set based on inspection of the data" (from
    #   their vignette section 2.1); I chose 20 because we know there are 20
    #   pathways a priori. We know nothing else about the data inside this loop.
    #   The proc.row and w.data parameters are set to yield MFA (see ?moa).
    geneMethyl_mgsa <- mogsa(
      x = multiOmics_ls,
      sup = geneSetMat_ls,
      nf = 20,
      proc.row = "center_ssq1",
      w.data = "lambda1"
    )
    # 1.34 seconds
    
    p.mat <- getmgsa(geneMethyl_mgsa, "p.val")
    mgsaCounts <- rowSums(p.mat < signif.level)
    # pathways x n matrix of p-values
    # sort(rowSums(p.mat < signif.level), decreasing = TRUE)
    # This yields the following count vector:
    # c(
    #   path7 = 41, path9 = 18, path15 = 18, path19 = 11, path6 = 9, path10 = 7,
    #   path11 = 6, path5 = 5, path16 = 4, path18 = 4, path20 = 4, path1 = 3,
    #   path3 = 3, path17 = 3, path14 = 2, path4 = 1, path8 = 1, path13 = 1,
    #   path2 = 0, path12 = 0
    # )
    # I don't know how to say which pathway is significant based on these
    #   counts. Perhaps a binomial probability with success 5%?
    
    ###  Binomial Probability  ###
    # E.g.: I want to find the probability that I flip 10 coins and get 2 or
    #   more heads
    # 1 - pbinom(2, size = 10, prob = 0.5)
    # Now I want to find the probability that, given that the null hypothesis
    #   is true, I find 41 or more "heads" for pathway 7
    pathwayPvals <- sapply(mgsaCounts, function(N){
      1 - pbinom(N, size = n.samples, prob = signif.level)
    })
    mogsaSignifPaths_df <- data.frame(
      terms = names(pathwayPvals),
      rawp = pathwayPvals,
      row.names = NULL
    )
    
    ###  Bootstrap the Background Counts  ###
    # bootPvals_mat <- sapply(1:10, function(i){
    #   
    #   pValsBoot_mat <- p.mat[
    #     , sample(seq_len(ncol(p.mat)), size = ncol(p.mat), replace = TRUE)
    #     ]
    #   mgsaCountsBoot <- rowSums(pValsBoot_mat < signif.level)
    #   binompVals <- sapply(mgsaCountsBoot, function(N){
    #     1 - pbinom(N, size = n.samples, prob = signif.level)
    #   })
    #   
    #   matrix(binompVals, nrow = 1)
    #   
    # })
    # rownames(bootPvals_mat) <- rownames(p.mat)
    # sort(rowMeans(bootPvals_mat))
    
    
    
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
    signifPathsG_logi <- sapply(
      seq_len(n.groups),
      isGroupDE,
      grpCard_int = geneGrpCard_int,
      DE_logi = exprGene_logi,
      pctDE = localParams_num[2] / 100
    )
    signifPathsG_int <- which(signifPathsG_logi)
    
    signifPathsG_df <- data.frame(
      terms = names(signifPathsG_logi),
      dePath = signifPathsG_logi
    )
    
    mogsaSignifPaths_df <- inner_join(mogsaSignifPaths_df, signifPathsG_df)
    
    rm(
      diffGenes_char, diffGenes_idx, exprGene_logi,
      signifPathsG_logi, signifPathsG_df
    )
    
    
    
    ######  Sensitivity and Specificity  ######################################
    
    pValues_num <- seq(0.001, 0.999, by = 0.005)
    
    geneSnS_df <- lapply(pValues_num, function(pVal){
      
      decisionTable(
        PathpVals_df = mogsaSignifPaths_df,
        truePaths = names(signifPathsG_int),
        type = "GeneExp",
        pValThresh = pVal
      )
      
    }) %>% 
      bind_rows
    
    multiOmicsSnS <-
      geneSnS_df %>% 
      arrange(pValue) %>% 
      mutate(TPR = TP / (TP + FN)) %>% 
      mutate(FPR = FP / (FP + TN)) %>% 
      mutate(ACC = (TP + TN) / n.groups) %>% 
      mutate(F1 = 2 * TP / (2 * TP + FP + FN)) %>% 
      mutate(
        MCC = (TP * TN - FP * FN) / 
          sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
      )
    
    # ggplot(multiOmicsSnS) +
    #   aes(x = pValue, y = MCC, group = DataType, colour = DataType) +
    #   geom_point(alpha = 0.5) +
    #   geom_vline(xintercept = 0.05) +
    #   geom_hline(yintercept = 0)
    # 
    # ggplot(multiOmicsSnS) +
    #   aes(x = pValue, y = ACC, group = DataType, colour = DataType) +
    #   geom_point(alpha = 0.5) +
    #   geom_vline(xintercept = 0.05)
    
    
    rm(geneSnS_df, pValues_num)
    
    
    
    ######  Simulation Experiment Output  #####################################
    run_int <- as.integer(str_replace_all(runDir, "[^0-9]", ""))
    multiOmicsSnS$run       <- run_int
    mogsaSignifPaths_df$run <- run_int
    
    out_ls <- list(
      fileNames = c(geneData_path, methylData_path),
      parameters = localParams_num,
      classification = multiOmicsSnS,
      multiOmicsSignif = mogsaSignifPaths_df
    )
    
    closeAllConnections()
    
    out_ls
    
    # END lapply() Runs
    
  })
  Sys.time() - a
  # 1.256394 min for first 15; 8.781544 min for 100 runs
  
  
  
  ######  Return  #############################################################
  confusion_df <-
    lapply(simResults_ls, `[[`, "classification") %>% 
    bind_rows() %>% 
    select(pValue, run, everything(), -DataType) %>% 
    arrange(pValue, run)
  mogsaPvalues_df <- 
    lapply(simResults_ls, `[[`, "multiOmicsSignif") %>% 
    bind_rows() %>% 
    select(run, everything())
  
  out_ls <- list(
    parameters = localParams_num,
    confusion = confusion_df,
    pValues = mogsaPvalues_df
  )
  
  outFile_char <- paste(
    "mogsa",
    paste(localParams_num, collapse = "_"),
    "all",
    sep = "_"
  )
  saveRDS(
    out_ls, file = paste0("results/sim_mogsa_20190612/", outFile_char, ".RDS")
  )
  
  # END for designs
  
}

Sys.time() - a

# For 12 design points iterating over 100 simulation runs:
# 16.61194 hours
