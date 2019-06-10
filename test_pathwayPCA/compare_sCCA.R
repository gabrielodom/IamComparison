# Test sCCA Multi-Omics
# Gabriel Odom
# 2019-06-10



######  Description and Overview  #############################################
# We have simulation results from sCCA under the same simulation parameters as
#   the pathwayPCA multi-omics sim. They are
#   - ngroups = 20
#   - p = 1600 (gene expression) & 2400 (DNA methylation)
#   - DE% = 10, 20, 50, 100%
#   - DEsize = 2x, 4x, 8x
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
# I'm not 100% sure on this. Should I be comparing the gene expression or the
#   DNA methylation? The counts for methylation will be out of 2400 instead of
#   1600. Also, n1, n2, n3, and n3 are based on the genes *identified as* DE by
#   the sCCA algorithm, not the truth. This should give me a p-value for each
#   pathway.
# EDIT: Lily said to use the 1600 for gene expression.

# Now that we have a p-value for each pathway, we also need to know the true
#   pathways. For this, we have to call back to the original data, and calculate
#   which pathways are truly DE. Then, finally, we can compare the decisions
#   made with the Fisher's Exact Test with the true pathways. I am expecting the
#   power of sCCA to be quite high (as Fisher's test can be sensitive to small
#   changes when you have very large counts), but I think the Type-I error of
#   sCCA might be inflated.


######  Setup  ################################################################
library(tidyverse)
library(pathwayPCA)


# The global parameters are
p.gene.exp <- 1600
p.methylation <- 2400
n.groups <- 20
geneGrpCard_int <- ceiling(p.gene.exp / n.groups)
methylGrpCard_int <- ceiling(p.methylation / n.groups)


# The design component of the file names are
designs_char <- c(
  "20gr_100p_80", "20gr_100p_40", "20gr_100p_20",
  "20gr_50p_80", "20gr_50p_40", "20gr_50p_20",
  "20gr_20p_80", "20gr_20p_40", "20gr_20p_20",
  "20gr_10p_80", "20gr_10p_40", "20gr_10p_20"
)
desg <- 1

localParams_num <- 
  str_split(designs_char[desg], "_") %>% 
  unlist() %>% 
  str_extract("[0-9]*") %>% 
  as.numeric()
names(localParams_num) <- c("numGroups", "pctDE", "pctEffectSize")

# Helper functions are
source("test_pathwayPCA/global_functions.R")



######  Pathway Collections  ##################################################
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


######  Directories  ##########################################################
# The data files all reside in
dataDir <- "results/GeneExp_Met_2DS/"

# Simulated data are in
simDir <- paste0(
  dataDir, "syntheticData/gr_20_perc_10_100_es_20_80_bkgrnd_0/RData/"
)

# sCCA Results are in
resDir <- paste0(dataDir, "synthetsCCA/supervised_exact-ES/files/")
runs <- paste0("run", 1:15, "/")

# a <- Sys.time()
res100_80_ls <- lapply(runs, function(runDir){

#



######  Inspect Current Comparison Results  ###################################

# run1_acc_arr <- readRDS(paste0(resDir, runDir, "sCCA_acc_arr.RData"))
# run1_auc_arr <- readRDS(paste0(resDir, runDir, "sCCA_auc_arr.RData"))
# run1_f1_arr <- readRDS(paste0(resDir, runDir, "sCCA_f1_arr.RData"))
# run1_FPR_arr <- readRDS(paste0(resDir, runDir, "sCCA_FPR_arr.RData"))
# run1_kappa_arr <- readRDS(paste0(resDir, runDir, "sCCA_kappa_arr.RData"))
# run1_mcc_arr <- readRDS(paste0(resDir, runDir, "sCCA_mcc_arr.RData"))
# run1_PW_mcc_arr <- readRDS(paste0(resDir, runDir, "sCCA_PW_mcc_arr.RData"))
# run1_TPR_arr <- readRDS(paste0(resDir, runDir, "sCCA_TPR_arr.RData"))
# 
# # These don't really have the information we need.
# rm(
#   run1_acc_arr, run1_auc_arr, run1_f1_arr, run1_FPR_arr, run1_kappa_arr,
#   run1_mcc_arr, run1_PW_mcc_arr, run1_TPR_arr
# )



######  Extract New Results  ##################################################
# We need to know which features sCCA claims are DE.

# sCCA Genes
geneRes_path <- paste0("sCCA_GeneExp_synth_", designs_char[desg], ".txt")
sCCAgenes_char <- scan(
  paste0(resDir, runDir, geneRes_path),
  what = character()
)

# sCCA Methylated Sites
methylRes_path <- paste0("sCCA_Methyl_synth_", designs_char[desg], ".txt")
sCCAmethyl_char <- scan(
  paste0(resDir, runDir, methylRes_path),
  what = character()
)



######  Fisher's Exact Test for All Pathways  #################################

# ###  Gene Expression  ###
# # Genes in pw1 & marked as DE
# n1 <- length(
#   intersect(
#     gene_PC$pathways$path1,
#     sCCAgenes_char
#   )
# )
# # Genes !(in pw1) & marked as DE
# n2 <- length(
#   intersect(
#     setdiff(allGenes_char, gene_PC$pathways$path1),
#     sCCAgenes_char
#   )
# )
# # Genes in pw1 & !(marked as DE)
# n3 <- length(
#   intersect(
#     gene_PC$pathways$path1,
#     setdiff(allGenes_char, sCCAgenes_char)
#   )
# )
# # Genes !(in pw1) & !(marked as DE)
# n4 <- length(
#   intersect(
#     setdiff(allGenes_char, gene_PC$pathways$path1),
#     setdiff(allGenes_char, sCCAgenes_char)
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
pathwayFisherExact(
  pathway = gene_PC$pathways$path1,
  allFeatures_char = allGenes_char,
  deFeatures_char = sCCAgenes_char
)


###  Apply Time  ###
sCCA_pVals <- sapply(
  gene_PC$pathways,
  pathwayFisherExact,
  allFeatures_char = allGenes_char,
  deFeatures_char = sCCAgenes_char
)



######  True DE Pathways  #####################################################
# Now we need to know which features are truly DE

###  Differentially-expressed genes  ###
geneDE_path <- paste0("GeneExp_synth_", designs_char[desg], "_Diff.RData")
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


# ###  Differentially-methylated sites  ###
# methylDE_path <- paste0("Methyl_synth_", designs_char[desg], "_Diff.RData")
# diffMethyl_char <- readRDS(paste0(simDir, runDir, methylDE_path))
# 
# diffMethyl_idx <- 
#   diffMethyl_char %>% 
#   str_replace("met_", "") %>% 
#   as.integer() %>% 
#   sort
# 
# methyl_logi <- seq_len(p.methylation) %in% diffMethyl_idx
# 
# # Which pathways were expressed?
# signifPathsM_int <- which(
#   sapply(
#     seq_len(n.groups),
#     isGroupDE,
#     grpCard_int = methylGrpCard_int,
#     DE_logi = methyl_logi,
#     pctDE = localParams_num[2] / 100
#   )
# )



######  Return  ###############################################################
data.frame(
  pVals = sCCA_pVals,
  dePath = seq_along(sCCA_pVals) %in% signifPathsG_int
)

})
# b <- Sys.time()
# b - a # 0.7188621 sec for 15 reps

names(res100_80_ls) <- paste0("run", 1:15)
res100_80_df <- bind_rows(res100_80_ls, .id = "run")
saveRDS(res100_80_df, file = "results/sim_sCCA_20190610/res_100_80.RDS")

res100_80_df %>% 
  mutate(claimDE = pVals < 0.05) %>% 
  mutate(T1err = claimDE & !dePath) %>% 
  mutate(T2err = (!claimDE) & dePath) %>% 
  group_by(run) %>% 
  summarise(
    pctT1 = mean(T1err),
    pctT2 = mean(T2err)
  )
