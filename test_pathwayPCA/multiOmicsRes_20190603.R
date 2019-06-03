# Simulation Run
# Gabriel Odom
# 2019-05-29

# This is a modified copy of test_pathwayPCA/multiOmicsRes_20190529.R, but this
#   version will analyse synthetic data with 0% background genes DE. The goal
#   is to find out how many of the pathways were correctly identified as
#   significant (measuring sensitivity and specificity), which pathways are
#   correctly identified as significant / not significant in *both* data sets, 
#   and the test size  / type 1 error for potential global tests (this last
#   metric is why we simulated an entirely new suite of data).


library(tidyverse)
library(pathwayPCA)



######  Global Functions  #####################################################
makeGroupIdx <- function(grp_idx, grpCard_int){
  
  pwy_idx <- seq.int(
    from = (grp_idx - 1) * grpCard_int + 1,
    to = grp_idx * grpCard_int
  )
  
}

makePath <- function(grp_idx, grpCard_int, pathType = c("ge_", "met_")){
  
  pathType <- match.arg(pathType)
  
  pwy_idx <- makeGroupIdx(grp_idx, grpCard_int)
  paste0(pathType, pwy_idx)
  
}

isGroupDE <- function(grp_idx, grpCard_int, DE_logi, pctDE){
  
  pwy_idx <- makeGroupIdx(grp_idx, grpCard_int)
  out_logi <- mean(DE_logi[pwy_idx]) >= pctDE
  names(out_logi) <- paste0("path", grp_idx)
  
  out_logi
  
}

decisionTable <- function(PathpVals_df, truePaths,
                          type, pValThresh = 0.01){
  
  # TP: marking a pathway as significant when it was designed as significant
  tpCount <- PathpVals_df %>% 
    filter(rawp <= pValThresh) %>% 
    select(terms) %>% 
    transmute(TP = terms %in% truePaths) %>%
    summarise(sum(TP)) %>%
    pull
  
  # FN: marking a pathway as not significant when it was designed as significant
  fnCount <- PathpVals_df %>% 
    filter(rawp > pValThresh) %>% 
    select(terms) %>% 
    transmute(TP = terms %in% truePaths) %>%
    summarise(sum(TP)) %>%
    pull
  
  # FP: marking a pathway as significant when it was designed not significant
  fpCount <- PathpVals_df %>% 
    filter(rawp <= pValThresh) %>% 
    select(terms) %>% 
    transmute(TP = !(terms %in% truePaths)) %>%
    summarise(sum(TP)) %>%
    pull
  
  # TN: marking a pathway as not significant when it was designed not significant
  tnCount <- PathpVals_df %>% 
    filter(rawp > pValThresh) %>% 
    select(terms) %>% 
    transmute(TP = !(terms %in% truePaths)) %>%
    summarise(sum(TP)) %>%
    pull
  
  data.frame(
    DataType = type,
    pValue = pValThresh,
    TP = tpCount,
    FN = fnCount,
    FP = fpCount,
    TN = tnCount,
    stringsAsFactors = FALSE
  )
  
}

pathSignif <- function(pathway, resp, omicsOut1, omicsOut2){
  
  path_df <- data.frame(
    Resp = resp,
    ome1 = getPathPCLs(omicsOut1, pathway)$PCs$V1,
    ome2 = getPathPCLs(omicsOut2, pathway)$PCs$V1
  )
  
  
  # Individual p-values
  ome1_mod <- glm(Resp ~ ome1, family = binomial, data = path_df)
  ome2_mod <- glm(Resp ~ ome2, family = binomial, data = path_df)
  
  pVals_df <- data.frame(
    ome1 = coef(summary(ome1_mod))[-1, 4],
    ome2 = coef(summary(ome2_mod))[-1, 4]
  )
  rownames(pVals_df) = "pValue"
  
  
  # global p-value
  null_mod <- glm(Resp ~ 1, family = binomial, data = path_df)
  path_mod <- glm(Resp ~ ., family = binomial, data = path_df)
  # glm(y ~ .) includes the intercept, so these models are nested.
  
  path_aov <- anova(null_mod, path_mod)
  LRpVal <- pchisq(path_aov[2, 4], df = path_aov[2, 3], lower.tail = FALSE)
  pVals_df$global <- LRpVal
  
  
  # Return
  pVals_df
  
}


######  Global Parameters  ####################################################
p.gene.exp <- 1600
p.methylation <- 2400
n.groups <- 20
geneGrpCard_int <- ceiling(p.gene.exp / n.groups)
methylGrpCard_int <- ceiling(p.methylation / n.groups)

whereDir <- 
  "results/GeneExp_Met_2DS/syntheticData/gr_20_perc_10_100_es_20_80_bkgrnd_0/RData/"
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
  "20gr_100p_80", "20gr_100p_40", "20gr_100p_20",
  "20gr_50p_80", "20gr_50p_40", "20gr_50p_20",
  "20gr_20p_80", "20gr_20p_40", "20gr_20p_20",
  "20gr_10p_80", "20gr_10p_40", "20gr_10p_20"
)
# desg <- 1
# designs_char[desg]
# # 20 groups; 100 percent of features DE in perturbed pathways; 80% effect size

a <- Sys.time()

for(desg in designs_char){

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



######  Runs Loop  ############################################################

# runDir <- runDirs_char[1]

for(runDir in runDirs_char){
  
  ######  Differentially-Expressed Genes  #######################################
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
  signifPathsG_int <- which(
    sapply(
      seq_len(n.groups),
      isGroupDE,
      grpCard_int = geneGrpCard_int,
      DE_logi = exprGene_logi,
      pctDE = localParams_num[2] / 100
    )
  )
  
  rm(diffGenes_char, diffGenes_idx, exprGene_logi)
  
  
  
  ######  Genes  ################################################################
  gene_mat <- readRDS(
    paste0(whereDir, runDir, geneData_path)
  )
  
  tGene_df <- data.frame(
    SampleIDs = colnames(gene_mat),
    Type = rep(c("tumour", "normal"), each = 100),
    t(gene_mat),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  
  ######  Gene OmicsCateg  ######################################################
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
  
  
  
  ######  Differentially-Methylated Genes  ######################################
  
  diffMethyl_char <- readRDS(
    paste0(whereDir, runDir, diffMethyl_path)
  )
  
  diffMethyl_idx <- 
    diffMethyl_char %>% 
    str_replace("met_", "") %>% 
    as.integer() %>% 
    sort
  
  methyl_logi <- seq_len(p.methylation) %in% diffMethyl_idx
  
  # Which pathways were expressed?
  signifPathsM_int <- which(
    sapply(
      seq_len(n.groups),
      isGroupDE,
      grpCard_int = methylGrpCard_int,
      DE_logi = methyl_logi,
      pctDE = localParams_num[2] / 100
    )
  )
  
  rm(diffMethyl_char, diffMethyl_idx, methyl_logi)
  
  
  
  ######  Methylation  ##########################################################
  methyl_mat <- readRDS(
    paste0(whereDir, runDir, methylData_path)
  )
  
  tMethyl_df <- data.frame(
    SampleIDs = colnames(methyl_mat),
    Type = rep(c("tumour", "normal"), each = 100),
    t(methyl_mat),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  
  ######  Methylation OmicsCateg  ###############################################
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
  
  
  
  ######  Sensitivity and Specificity  ##########################################
  
  pValues_num <- seq(0.001, 0.999, by = 0.005)
  
  geneSnS_df <- lapply(pValues_num, function(pVal){
    
    decisionTable(
      PathpVals_df = getPathpVals(gene_aespcOut),
      truePaths = names(signifPathsG_int),
      type = "GeneExp",
      pValThresh = pVal
    )
    
  }) %>% 
    bind_rows
  
  methylSnS_df <- lapply(pValues_num, function(pVal){
    
    decisionTable(
      PathpVals_df = getPathpVals(methyl_aespcOut),
      truePaths = names(signifPathsM_int),
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
  
  
  rm(geneSnS_df, methylSnS_df, pValues_num)
  
  
  
  ######  Global p-Value  #######################################################
  
  multiOmicspVals_df <- lapply(
    getPathwayCollection(gene_Omics)$TERMS,
    pathSignif,
    resp = getResponse(gene_Omics),
    omicsOut1 = gene_aespcOut,
    omicsOut2 = methyl_aespcOut
  ) %>% 
    bind_rows(.id = "pathway") %>% 
    mutate(DEpath = pathway %in% names(signifPathsG_int)) %>% 
    arrange(global)
  
  
  
  ######  Simulation Experiment Output  #########################################
  simResults_ls <- list(
    fileNames = c(geneData_path, methylData_path),
    parameters = localParams_num,
    classification = multiOmicsSnS,
    multiOmicsSignif = multiOmicspVals_df,
    geneOmics = gene_Omics,
    geneResults = gene_aespcOut,
    geneTruth = signifPathsG_int,
    methylOmics = methyl_Omics,
    methylResults = methyl_aespcOut,
    methylTruth = signifPathsG_int
  )
  
  out_path <- "results/sim_20190603/"
  outFile_char <- paste(
    "multiOmics",
    paste(localParams_num, collapse = "_"),
    str_extract(runDir, "(\\d)+"), 
    sep = "_"
  )
  
  saveRDS(simResults_ls, file = paste0(out_path, outFile_char, ".RDS"))
  
  # END for() Runs
  
}

# END for designs

}

Sys.time() - a

# For iterating over 100 simulation runs (in hours):
# 1.402108, 1.398272, 1.411486, 1.401245, 1.399975, 1.412866, 1.404951, 1.41051,
#   1.411761


# For iterating over 9 design points:
# ???, 7.646102, 7.734493, 7.804814, 7.648272, 7.674367, 7.408685, 7.87128,
#   7.945855