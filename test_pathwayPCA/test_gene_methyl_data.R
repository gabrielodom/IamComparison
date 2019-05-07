library(tidyverse)
# devtools::install_github("gabrielodom/pathwayPCA", ref = "stable_3_5")
library(pathwayPCA)


###  Pathways  ###
diffMethyl_char <- readRDS(
  "results/GeneExp_Met_2DS/syntheticData/gr_20_perc_10_100_es_20_80/RData/run1/Methyl_synth_20gr_100p_80_Diff.RData"
)
sort(diffMethyl_char)
diffMethyl_idx <- sort(
  as.integer(
    str_replace(diffMethyl_char, "met_", "")
  )
)
plot(diffMethyl_idx)
diffMethyl_idx[1:121]
diffMethyl_idx[2:122]
# There are 80 genes in this pathway: 2400 / 20 = 120. Therefore, all of these
#   pathways are the same size. The parameter .percentage.pw controls the 
#   percentage of genes in each pathway are artificially differentially-
#   expressed (for this run, we see "...perc_10_100...", so the chosen pathways
#   are 100% expressed).
percentage.pw <- 1
exprMethyl_logi <- 1:2400 %in% diffMethyl_idx
plot(exprMethyl_logi)

# Which pathways were expressed?
signifPathsM_int <- which(
  sapply(1:20, function(p){
    
    pLen_int <- 120
    pwy_int <- seq.int(from = (p - 1) * pLen_int + 1, p * pLen_int)
    
    out_num <- mean(exprMethyl_logi[pwy_int]) >= percentage.pw
    names(out_num) <- paste0("path", p)
    
    out_num
    
  })
)


methylPaths_ls <- lapply(1:20, function(p){
  
  pLen_int <- 120
  pwy_int <- seq.int(from = (p - 1) * pLen_int + 1, p * pLen_int)
  paste0("met_", pwy_int)
  
})
names(methylPaths_ls) <- paste0("path", 1:length(methylPaths_ls))

methyl_PC <- CreatePathwayCollection(
  sets_ls = methylPaths_ls,
  TERMS = names(methylPaths_ls)
)



###  Methyl  ###
methyl_mat <- readRDS(
  "results/GeneExp_Met_2DS/syntheticData/gr_20_perc_10_100_es_20_80/RData/run1/Methyl_synth_20gr_100p_80.RData"
)

dim(methyl_mat)
colnames(methyl_mat)

tMethyl_df <- data.frame(
  SampleIDs = colnames(methyl_mat),
  Type = rep(c("tumour", "normal"), each = 100),
  t(methyl_mat),
  stringsAsFactors = FALSE
)


###  Test OmicsCateg  ###
methyl_Omics <- CreateOmics(
  assayData_df = tMethyl_df[, -2],
  pathwayCollection_ls = methyl_PC,
  response = tMethyl_df[, 1:2],
  respType = "categ"
)

methyl_aespcOut <- AESPCA_pVals(
  methyl_Omics,
  parallel = TRUE,
  numCores = 4,
  adjustment = "BH"
)

getPathpVals(methyl_aespcOut)
signifPathsM_int
plot(exprMethyl_logi)
# We find exactly the correct pathways, but we also find pathway 19.