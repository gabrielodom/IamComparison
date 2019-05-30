# Simulation Run
# Gabriel Odom
# 2019-05-29

# This is a slimmed-down version of test_pathwayPCA/test_gene_expr_data.R and 
#   test_pathwayPCA/test_gene_methyl_data.R. The goal is to find out how many
#   of the pathways were correctly identified as significant (measuring
#   sensitivity and specificity) and which pathways are correctly identified as
#   significant / not significant in *both* data sets.


library(tidyverse)
library(pathwayPCA)

whichRun_int <- 1
where_char <- paste0(
  "results/GeneExp_Met_2DS/syntheticData/gr_20_perc_10_100_es_20_80/RData/run",
  whichRun_int,
  "/"
)
geneDataFiles_char <- list.files(where_char, pattern = "^GeneExp.*Diff.RData$")
methylDataFiles_char <- list.files(where_char, pattern = "^Methyl.*Diff.RData$")


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
                          type, params, FDR_num = 0.05){
  
  # TP: marking a pathway as significant when it was designed as significant
  tpCount <- PathpVals_df %>% 
    filter(FDR_BH <= FDR_num) %>% 
    select(terms) %>% 
    transmute(TP = terms %in% truePaths) %>%
    summarise(sum(TP)) %>%
    pull
  
  # FN: marking a pathway as not significant when it was designed as significant
  fnCount <- PathpVals_df %>% 
    filter(FDR_BH > FDR_num) %>% 
    select(terms) %>% 
    transmute(TP = terms %in% truePaths) %>%
    summarise(sum(TP)) %>%
    pull
  
  # FP: marking a pathway as significant when it was designed not significant
  fpCount <- PathpVals_df %>% 
    filter(FDR_BH <= FDR_num) %>% 
    select(terms) %>% 
    transmute(TP = !(terms %in% truePaths)) %>%
    summarise(sum(TP)) %>%
    pull
  
  # TN: marking a pathway as not significant when it was designed not significant
  tnCount <- PathpVals_df %>% 
    filter(FDR_BH > FDR_num) %>% 
    select(terms) %>% 
    transmute(TP = !(terms %in% truePaths)) %>%
    summarise(sum(TP)) %>%
    pull
  
  data.frame(
    DataType = type,
    data.frame(t(params)),
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
  
  null_mod <- glm(Resp ~ 1, family = binomial, data = path_df)
  path_mod <- glm(Resp ~ ., family = binomial, data = path_df)
  
  path_aov <- anova(null_mod, path_mod)
  LRpVal <- pchisq(path_aov[2, 4], df = path_aov[2, 3], lower.tail = FALSE)
  
  pVals_mat <- t(coef(summary(path_mod)))[4, -1, drop = FALSE]
  pVals_df <- as.data.frame(pVals_mat, row.names = "pValue")
  pVals_df$global <- LRpVal
  
  pVals_df
  
}


######  Global Parameters  ####################################################
p.gene.exp <- 1600
p.methylation <- 2400
n.groups <- 20
geneGrpCard_int <- ceiling(p.gene.exp / n.groups)
methylGrpCard_int <- ceiling(p.methylation / n.groups)


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



######  Differentially-Expressed Genes  #######################################
diffGenes_path <- geneDataFiles_char[3]

# 20 groups; 100 percent of features DE in perturbed pathways; 80% effect size
localParams_num <- 
  str_split(diffGenes_path, "_") %>% 
  unlist() %>% 
  `[`(3:5) %>% 
  str_extract("[0-9]*") %>% 
  as.numeric()
names(localParams_num) <- c("numGroups", "pctDE", "pctEffectSize")


diffGenes_char <- readRDS(
  paste0(where_char, diffGenes_path)
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
geneData_path <- str_remove(geneDataFiles_char[3], "_Diff")
gene_mat <- readRDS(
  paste0(where_char, geneData_path)
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

rm(gene_mat, tGene_df, gene_PC)



######  Differentially-Methylated Genes  ######################################
diffMethyl_path <- methylDataFiles_char[3]

diffMethyl_char <- readRDS(
  paste0(where_char, diffMethyl_path)
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
methylData_path <- str_remove(methylDataFiles_char[3], "_Diff")
methyl_mat <- readRDS(
  paste0(where_char, methylData_path)
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

rm(methyl_mat, tMethyl_df, methyl_PC)



######  Sensitivity and Specificity  ##########################################
geneSnS_df <- decisionTable(
  PathpVals_df = getPathpVals(gene_aespcOut),
  truePaths = names(signifPathsG_int),
  type = "GeneExp",
  params = localParams_num
)

methylSnS_df <- decisionTable(
  PathpVals_df = getPathpVals(methyl_aespcOut),
  truePaths = names(signifPathsM_int),
  type = "Methylation",
  params = localParams_num
)

multiOmicsSnS <- rbind(geneSnS_df, methylSnS_df)

rm(geneSnS_df, methylSnS_df)



######  Global p-Value  #######################################################

multiOmicspVals_df <- lapply(
  getPathwayCollection(gene_Omics)$TERMS,
  pathSignif,
  resp = getResponse(gene_Omics),
  omicsOut1 = gene_aespcOut,
  omicsOut2 = methyl_aespcOut
) %>% 
  bind_rows(.id = "pathway")



######  Simulation Experiment Output  #########################################
simResults_ls <- list(
  fileNames = c(geneData_path, methylData_path),
  classification = multiOmicsSnS,
  geneOmics = gene_Omics,
  geneResults = gene_aespcOut,
  geneTruth = signifPathsG_int,
  methylOmics = methyl_Omics,
  methylResults = methyl_aespcOut,
  methylTruth = signifPathsG_int
)

out_path <- "results/sim_20190529/"
outFile_char <- paste(
  "multiOmics",
  paste(localParams_num, collapse = "_"),
  whichRun_int,
  sep = "_"
)

saveRDS(simResults_ls, file = paste0(out_path, outFile_char, ".RDS"))
