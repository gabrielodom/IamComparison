resDir <- "results/GeneExp_Met_2DS/synthetsCCA/supervised_exact-ES/RData/"
runDir <- "run1/"
run1_acc_arr <- readRDS(paste0(resDir, runDir, "sCCA_acc_arr.RData"))
run1_auc_arr <- readRDS(paste0(resDir, runDir, "sCCA_auc_arr.RData"))
run1_f1_arr <- readRDS(paste0(resDir, runDir, "sCCA_f1_arr.RData"))
run1_FPR_arr <- readRDS(paste0(resDir, runDir, "sCCA_FPR_arr.RData"))
run1_kappa_arr <- readRDS(paste0(resDir, runDir, "sCCA_kappa_arr.RData"))
run1_mcc_arr <- readRDS(paste0(resDir, runDir, "sCCA_mcc_arr.RData"))
run1_PW_mcc_arr <- readRDS(paste0(resDir, runDir, "sCCA_PW_mcc_arr.RData"))
run1_TPR_arr <- readRDS(paste0(resDir, runDir, "sCCA_TPR_arr.RData"))

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
# Now that we have a p-value for each pathway, we also need to know the true
#   pathways. For this, we have to call back to the original data, and calculate
#   which pathways are truly DE. Then, finally, we can compare the decisions
#   made with the Fisher's Exact Test with the true pathways. I am expecting the
#   power of sCCA to be quite high (as Fisher's test can be sensitive to small
#   changes when you have very large counts), but I think the Type-I error of
#   sCCA might be inflated.