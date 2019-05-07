library(tidyverse)
devtools::install_github("gabrielodom/pathwayPCA", ref = "stable_3_5")
library(pathwayPCA)


###  Pathways  ###
diffGenes_char <- readRDS(
  "results/GeneExp_Met_2DS/syntheticData/gr_20_perc_10_100_es_20_80/RData/run1/GeneExp_synth_20gr_100p_80_Diff.RData"
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



###  Genes  ###
gene_mat <- readRDS(
  "results/GeneExp_Met_2DS/syntheticData/gr_20_perc_10_100_es_20_80/RData/run1/GeneExp_synth_20gr_100p_80.RData"
)

dim(gene_mat)
colnames(gene_mat)

tGene_df <- data.frame(
  SampleIDs = colnames(gene_mat),
  Type = rep(c("tumour", "normal"), each = 100),
  t(gene_mat),
  stringsAsFactors = FALSE
)


###  Test OmicsCateg  ###
gene_Omics <- CreateOmics(
  assayData_df = tGene_df[, -2],
  pathwayCollection_ls = gene_PC,
  response = tGene_df[, 1:2],
  respType = "categ"
)

gene_aespcOut <- AESPCA_pVals(
  gene_Omics,
  parallel = TRUE,
  numCores = 4,
  adjustment = "BH"
)

getPathpVals(gene_aespcOut)
signifPathsG_int
# We find exactly the correct pathways, and only the correct pathways.