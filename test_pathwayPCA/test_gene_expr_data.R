library(tidyverse)
# devtools::install_github("gabrielodom/pathwayPCA", ref = "stable_3_5")
# BiocManager::install("pathwayPCA")
library(pathwayPCA)

where_char <-
  "results/GeneExp_Met_2DS/syntheticData/gr_20_perc_10_100_es_20_80/RData/run"
whichRun_int <- 1


######  Pathways  #############################################################
diffGenes_path <- "GeneExp_synth_20gr_100p_80_Diff.RData"
# 20 groups; 100 percent of features DE in perturbed pathways; 80% effect size

diffGenes_char <- readRDS(
  paste0(where_char, whichRun_int, "/", diffGenes_path)
)

sort(diffGenes_char)
diffGenes_idx <- sort(
  as.integer(
    str_replace(diffGenes_char, "ge_", "")
  )
)
plot(diffGenes_idx)
diffGenes_idx[1:81]
# There are 80 genes in this pathway: 1600 / 20 = 80. Therefore, all of these
#   pathways are the same size. The parameter .percentage.pw controls the 
#   percentage of genes in each pathway are artificially differentially-
#   expressed (for this run, we see "...perc_10_100...", so the chosen pathways
#   are 100% expressed).
percentage.pw <- 1
exprGene_logi <- 1:1600 %in% diffGenes_idx
plot(exprGene_logi)

# Which pathways were expressed?
signifPathsG_int <- which(
  sapply(1:20, function(p){
    
    pLen_int <- 80
    pwy_int <- seq.int(from = (p - 1) * pLen_int + 1, p * pLen_int)
    
    out_num <- mean(exprGene_logi[pwy_int]) >= percentage.pw
    names(out_num) <- paste0("path", p)
    
    out_num
    
  })
)


genePaths_ls <- lapply(1:20, function(p){
  
  pLen_int <- 80
  pwy_int <- seq.int(from = (p - 1) * pLen_int + 1, p * pLen_int)
  paste0("ge_", pwy_int)
  
})
names(genePaths_ls) <- paste0("path", 1:length(genePaths_ls))

gene_PC <- CreatePathwayCollection(
  sets_ls = genePaths_ls,
  TERMS = names(genePaths_ls)
)



######  Genes  ################################################################
geneData_path <- "GeneExp_synth_20gr_100p_80.RData"
gene_mat <- readRDS(
  paste0(where_char, whichRun_int, "/", geneData_path)
)

dim(gene_mat)
colnames(gene_mat)

tGene_df <- data.frame(
  SampleIDs = colnames(gene_mat),
  Type = rep(c("tumour", "normal"), each = 100),
  t(gene_mat),
  row.names = NULL,
  stringsAsFactors = FALSE
)


######  OmicsCateg  ###########################################################
gene_Omics <- CreateOmics(
  assayData_df = tGene_df[, -2],
  pathwayCollection_ls = gene_PC,
  response = tGene_df[, 1:2],
  respType = "categ"
)

system.time(
  gene_aespcOut <- AESPCA_pVals(
    gene_Omics,
    parallel = TRUE,
    numCores = 5,
    adjustment = "BH"
  )
)
# Time in seconds
# 29.47 for 2 cores; 22.23 for 4; 21.62 for 5; 24.28 for 6; 26.48 for 8

pVals_df <- getPathpVals(gene_aespcOut)
pVals_df
signifPathsG_int
# We find exactly the correct pathways, and only the correct pathways.



######  Sensitivity and Specificity  ##########################################

# TP: marking a pathway as significant when it was designed as significant
tpCount <- pVals_df %>% 
  filter(FDR_BH <= 0.05) %>% 
  select(terms) %>% 
  transmute(TP = terms %in% names(signifPathsG_int)) %>%
  summarise(sum(TP)) %>%
  pull
  
# FN: marking a pathway as not significant when it was designed as significant
fnCount <- pVals_df %>% 
  filter(FDR_BH > 0.05) %>% 
  select(terms) %>% 
  transmute(TP = terms %in% names(signifPathsG_int)) %>%
  summarise(sum(TP)) %>%
  pull

# FP: marking a pathway as significant when it was designed not significant
fpCount <- pVals_df %>% 
  filter(FDR_BH <= 0.05) %>% 
  select(terms) %>% 
  transmute(TP = !(terms %in% names(signifPathsG_int))) %>%
  summarise(sum(TP)) %>%
  pull

# TN: marking a pathway as not significant when it was designed not significant
tnCount <- pVals_df %>% 
  filter(FDR_BH > 0.05) %>% 
  select(terms) %>% 
  transmute(TP = !(terms %in% names(signifPathsG_int))) %>%
  summarise(sum(TP)) %>%
  pull

sns_df <- data.frame(TP = tpCount, FN = fnCount, FP = fpCount, TN = tnCount)
sns_df


######  Simulation Experiment Output  #########################################
params_int <-
  str_split(geneData_path, "_") %>%
  unlist %>% 
  `[`(3:5) %>% 
  str_extract("[0-9]*") %>% 
  as.integer()
names(params_int) <- c("eqGroups", "pct_DE", "%effectSize")

geneSimResults_ls <- list(
  simParams = params_int,
  fileNames = c(diffGenes_path, geneData_path),
  Omics = gene_Omics,
  Results = gene_aespcOut,
  truePaths = signifPathsG_int,
  classification = sns_df
)

out_path <- "results/sim_results_20190529/"
outFile_char <- paste(
  "genes",
  paste(params_int, collapse = "_"),
  whichRun_int,
  sep = "_"
)

saveRDS(geneSimResults_ls, file = paste0(out_path, outFile_char, ".RDS"))
