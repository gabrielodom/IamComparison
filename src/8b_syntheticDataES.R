#######################################################################
##### file: 8b_syntheticDataES.R                                  #####
##### input: parameters defined in 8a_synthetParameter.R          #####
##### output: synthetic datasets (as *.txt and *.RData) with      #####
#####   abberation level based on relative detectable effect      #####
##### packages:                                                   #####
##### author: B. Pucher                                           #####
##### date created: 04/10/2017                                    #####
##### last change: 09/11/2017                                     #####
#######################################################################
rm(list=ls())
ls.str(all.names = TRUE)

#install.packages("gplots")
library("gplots")

#######################################################################
#####                                                             #####
##### FUNCTIONS                                                   #####   
#####                                                             #####
#######################################################################

#####------------------------------------------------------------------
# simulate gene expression data from Gamma-distribution Wright and Simon (2003)
#####------------------------------------------------------------------
generateExpLevelsGamma = function(n_row, n_col){
  
  # draw gene expression values from normal distribution with 
  # sd of a gene drawn from inverse Gamma distribution 
  sigma_sq = 1/rgamma(n_row, shape = 3, scale = 1)
  ge_mtx = matrix(rnorm(n_row*n_col, mean = 0, sd = sqrt(sigma_sq)),
                  nrow = n_row, byrow = FALSE)
  
  return(ge_mtx) 
}

#####------------------------------------------------------------------
# Induce differential expression of features in diff_idx rows in 
# the dataset by adding a constant value c_de. Noise is also added.  
#####------------------------------------------------------------------
induceDE = function(basic_mtx, diff_idx, c, noise = 0.1, type){
  
  de_mtx = basic_mtx
  up_down = sample(c(1,-1), length(diff_idx), replace = TRUE)
  # difference of DE features between groups 
  # varies between features but is the same within a feature (for all samples)
  c_diff = rnorm(length(diff_idx), c, noise)

  de_up_down = matrix(up_down * c_diff, nrow = length(diff_idx), ncol = n_tu,
                                byrow = FALSE)
  
  de_mtx[diff_idx, 1:n_tu] =
    basic_mtx[diff_idx, 1:n_tu] + de_up_down
  colnames(de_mtx) = c(paste0(seq(1:n_tu), "_tu"), 
                       paste0(seq(1:n_no), "_no"))
  rownames(de_mtx) = paste0(type, "_", seq(1:nrow(de_mtx)))
  
  positives = rownames(de_mtx)[diff_idx]
  return(list(DE_mtx = de_mtx, Diff = positives))
  
} 

#####------------------------------------------------------------------
# simulate independent methylation levels (Li et al. 2015)
#####------------------------------------------------------------------
generateBetaMatrix = function(n_row, n_col){
  
  beta_mtx = matrix(sample(c(rbeta(0.1*n_row*n_col, 0.5, 5),
                             rbeta(0.9*n_row*n_col, 5, 0.5))), nrow = n_row)

  return(beta_mtx)
}

#######################################################################
#####                                                             #####
##### MAIN SECTION                                                #####
#####                                                             #####
#######################################################################
current.run = .current.synthet
sub.dir.name = file.path("syntheticData", current.run)
source(file.path(.src.dir, "setSubDirPath.R"))

log.con = file(file.path(sub.dir.files, "syntheticDataLog.txt"))
sink(file = log.con, type = "output")
flush(log.con)
cat(paste("---------------", Sys.time(), "---------------\n"))
cat("Start synthetic data generation...\n")

n_tu = .n.tumor; n_no = .n.normal # number of samples
n = n_tu+n_no
cat("Total number of samples for simulation: ", n, "\n")

# gene expression/RNA level dataset
m_ge = .m.gene.exp 
cat("Number of features in the gene expression dataset: ", m_ge, "\n")

# DNA methylation dataset
m_met = .m.methylation # number of features
cat("Number of features in the DNA methylation dataset: ", m_met, "\n")

# number of groups (pathways)
n_groups = .n.groups
cat("Number of pathways in Datasets: ", n_groups, "\n")

# number of perturbed pathways
n_groups_d = round(n_groups*.p.groups.d)
cat("Number of perturbed pathways in Datasets: ", n_groups_d, "\n")

# proportion of differential profiles in perturbed pathways
p_de_pw = .percentage.pw
p_dm_pw = .percentage.pw

# proportion of additional random differential profiles
p_random = .p.random

# proportion of differential profiles in total
p_de = 1/.n.groups*(.n.groups*.p.groups.d)*.percentage.pw
p_dm = 1/.n.groups*(.n.groups*.p.groups.d)*.percentage.pw
 
# relative effect to add for differential expression 
pt_de = .effect.size
pt_dm = .effect.size

# standard deviation (across features) of relative effect added to DE features
c_de_sd = 0.1
c_dm_sd = 0.1

# power to detect effect
pow = 0.8
alpha = 0.05

cat("Percentage of differential profiles (accumulation from perturbed pathway):
   gene expression ", 1/.n.groups*(.n.groups*.p.groups.d)*.percentage.pw, "%
   DNA methylation ", 1/.n.groups*(.n.groups*.p.groups.d)*.percentage.pw, "%\n")

cat("Percentage of additional randomly spread differential profiles:", p_random*100, "%\n")

cat("Relative effect used for simulation of DE features:", .effect.size, "\n")
cat("Standard deviation of relative effect across feautures: GE: ", c_de_sd, " MET: ", c_dm_sd, "\n")

runs = .stats.runs
cat("Create", runs, "datasets with each parameter set to achive valid",
    "statistics of AUC.\n")

sub.dir.figures.tmp = sub.dir.figures
sub.dir.files.tmp = sub.dir.files
sub.dir.RData.tmp = sub.dir.RData

# Seeds used in each stats run
stats.seed = .stats.seed

cat("Start data simulation ... \n")

# create example plots in first stats run
create.example.plot = TRUE

#####------------------------------------------------------------------
# simulation of gene expression data (Oana)
#####------------------------------------------------------------------
cat("Simulate gene expression dataset...  \n")

for(run in 1:runs){
  
  # reset sub.dir.paths
  sub.dir.figures = sub.dir.figures.tmp 
  sub.dir.files = sub.dir.files.tmp
  sub.dir.RData = sub.dir.RData.tmp
  
  # set sub.dir.paths to run dir
  sub.dir.figures = file.path(sub.dir.figures, paste0("run", run))
  sub.dir.files = file.path(sub.dir.files, paste0("run", run))
  sub.dir.RData = file.path(sub.dir.RData, paste0("run", run))  
  
  # create new directory for run
  dir.create(sub.dir.figures, showWarnings = FALSE)
  dir.create(sub.dir.files, showWarnings = FALSE)
  dir.create(sub.dir.RData, showWarnings = FALSE)
  
  # determine which pathways are perturbed 
  set.seed(stats.seed[run])
  groups_d = sample(n_groups, n_groups_d, replace = FALSE)
  cat("Stats run: ", run, " Pathways simulated as pertubed:", groups_d, "\n")
  # report seed used for current stats run
  cat("Seed used:", stats.seed[run], "\n")

  # generate basic gene expression matrix
  data_ge = generateExpLevelsGamma(n_row = m_ge, n_col = n)
  # plot example heatmap
  if(run == 1 & create.example.plot == TRUE) {
    png(file.path(sub.dir.figures, "GeneExp_synth_basic.png"), width = 1000, height = 1000)
    heatmap.2(data_ge, dendrogram = "none", trace = "none", Rowv = F,
              Colv = F, labRow = "", main = "Basic Gene Expression Dataset")
    invisible(dev.off())
  }
  
  # do not add additional noise to the data matrix
  data_ge_noise = data_ge 

  GE_PW_idx = factor(paste0("gr_", sort(rep(c(1:n_groups), m_ge/n_groups))))
  saveRDS(GE_PW_idx, file.path(sub.dir.RData, "GE_PWs.RData"))
  
  # Calculate detectable effect in each feature.
  # Effect (diff in means) detectable in t-test of two independent normal distibuted samples 
  # with approximately equal variances. (Sachs p. 385 or Bland p. 336 or Held p. 366) 
  # detectable effect delta
  effect_ge = apply(data_ge_noise, MARGIN = 1, function(x) {
    power.t.test(delta = NULL, n = (2*n_tu*n_no)/(n_tu+n_no), sig.level = alpha,
                 power = pow, sd = sd(x), type = "two.sample", alternative = "two.sided")$delta
  })

  effect_size_ge = effect_ge/apply(data_ge_noise, MARGIN = 1, sd)
  cat("Realtive detectable effect in each GE feature:", median(effect_size_ge), "\n")

  # plot histogram of detectable effect 
  png(file.path(sub.dir.figures, "Effect_size_GeneExp.png"))
  hist(effect_ge, main = paste("Effect detectable with n =", n), 
       breaks = 50, xlab = "Difference of Mean")
  invisible(dev.off())
  
  # determine factor for values to add to induce DE
  c_de = .effect.size
  dx = c_de/median(effect_size_ge)
  c_de_abs = dx*median(effect_ge)
  cat("Values to add for a median effect of", c_de, ":", round(c_de_abs, 3), "\n")
  
  # determine absolut standard deviation of effect 
  dx_sd = c_de_sd/median(effect_size_ge)
  c_de_sd_abs = dx_sd*median(effect_ge)

  
  if(run == 1 & create.example.plot == TRUE){
    # plot empirical distribution function of detectable effect size
    # and indicate percentiles of effect to add for DE
    png(file.path(sub.dir.figures, "GeneExp_Fn_effect.png"), width = 2000, height = 1000, res = 200)
    plot(ecdf(effect_ge), main = paste("Cumulative distribution of detectable Effect size (Power:", pow, ")"),
         xlab = "Difference of Mean", xaxt = "n")
    axis(side = 1, at = c(round(min(effect_ge),1), round(c_de_abs, 3), round(max(effect_ge), 1)),
         tick = c(round(min(effect_ge),1), round(c_de_abs, 3), max(round(max(effect_ge), 1), max(c_de))))
    abline(v = c_de_abs)
    text(x = c_de_abs, y = 0.8, labels = c_de, pos = 4)
    invisible(dev.off())
  }

  for(k in 1:length(n_groups)){
    
    # create loading matrix representing gene-groups/pathways
    if(m_ge %% n_groups[k] == 0 & m_met %% n_groups[k] == 0){

      # create datasets with different proportions of differentially expressed features
      for(i in 1: length(p_de_pw)){
        for(j in 1:length(c_de)){
          
          pw_de_idx = c(sapply(groups_d, function(x) {
            # within perturbed pathways
            sample(which(GE_PW_idx %in% paste0("gr_", x)), size = p_de_pw[i]/100*m_ge/n_groups, 
                   replace = FALSE)
          }))
          
          # randomly spread
          ran_de_idx = sample(1:m_ge, size = round(m_ge*p_random), replace = FALSE)
          
          diff_idx = unique(c(pw_de_idx, ran_de_idx))
            
          GeneExp_DE = induceDE(data_ge_noise, diff_idx, c = c_de_abs[j],
                                noise = c_de_sd_abs, type = "ge") 
          data_ge_de = GeneExp_DE$DE_mtx
          
          fname = paste("GeneExp_synth_", n_groups[k], "gr_", p_de_pw[i], "p_", pt_de[j]*100, sep = "")
          saveRDS(data_ge_de, file = file.path(sub.dir.RData, paste(fname, "RData", sep = ".")))
          #write.table(data_ge_de, file.path(sub.dir.files, paste(fname, "txt", sep = ".")), quote = FALSE)
          if(run == 1 & create.example.plot == TRUE) {
            png(file.path(sub.dir.figures, paste(fname, "png", sep = ".")), width = 1000, height = 1000)
            heatmap.2(data_ge_de, dendrogram = "none", trace = "none", Rowv = T, Colv = F,
                      labRow = "", main = paste("\n\nDE Genes in perturbed pathways:", p_de_pw[i],
                                                "% \n DE Genes in total: ", p_de[i],
                                                "% \nMean difference between groups in DE Genes:",
                                                "\npercentile: ", pt_de[j]*100, "/ value: ", c_de[j]))
                     invisible(dev.off())
          }
          
          positives = GeneExp_DE$Diff
          saveRDS(positives, file = file.path(sub.dir.RData, paste0(fname, "_Diff", ".RData")))
          write.table(positives, file.path(sub.dir.files, paste0(fname, "_Diff", ".txt")), 
                      row.names = FALSE, quote = FALSE)
          
        }
      }
    } else cat("Not all Features assigned to a Group! Quit.")
  }

  # assess whether DE can be detected with t-test
  ttest_result_ge = list()
  
  for(k in 1:length(n_groups)){
    for(i in 1: length(p_de_pw)){
      for(j in 1:length(pt_de)){
        
        fname = paste("GeneExp_synth_", n_groups[k], "gr_", p_de_pw[i], "p_", pt_de[j]*100, sep = "")
        data_ge_de = readRDS(file.path(sub.dir.RData, paste(fname, "RData", sep = ".")))
        positives = readRDS(file.path(sub.dir.RData, paste0(fname, "_Diff", ".RData")))
          
        ttest_p = apply(data_ge_de, 1, function(x) {
          t.test(x[1:n_tu], x[(n_tu+1):n], var.equal = TRUE, paired = FALSE,
                          alternative = "two.sided")$p.value })
        ttest_pFDR = p.adjust(ttest_p, method = "BH")
        positives_t = list(names(which(ttest_pFDR <= alpha)))
        names(positives_t) = fname
        ttest_result_ge = c(ttest_result_ge, positives_t)
      }
    }
  }
  
  saveRDS(ttest_result_ge, file = file.path(sub.dir.RData, "ttest_result_GeneExp.RData"))   
}


#####------------------------------------------------------------------
# Simulation of Methylation Data
#####------------------------------------------------------------------    
cat("Simulate DNA-Methylation dataset... \n")

for(run in 1:runs){
  
  # reset sub.dir.paths
  sub.dir.figures = sub.dir.figures.tmp 
  sub.dir.files = sub.dir.files.tmp
  sub.dir.RData = sub.dir.RData.tmp
  
  # set sub.dir.paths to run dir
  sub.dir.figures = file.path(sub.dir.figures, paste0("run", run))
  sub.dir.files   = file.path(sub.dir.files, paste0("run", run))
  sub.dir.RData   = file.path(sub.dir.RData, paste0("run", run))  
  
  # create new directory for run
  dir.create(sub.dir.figures, showWarnings = FALSE)
  dir.create(sub.dir.files, showWarnings = FALSE)
  dir.create(sub.dir.RData, showWarnings = FALSE)
  
  # determine which pathways are perturbed (same as GE data)
  set.seed(stats.seed[run])
  groups_d = sample(n_groups, n_groups_d, replace = FALSE)
  cat("Stats run: ", run, " Pathways simulated as pertubed:", groups_d, "\n")
  cat("Seed used:", stats.seed[run], "\n")
  
  # draw independent methylation levels from mixed Beta-distribution 
  # (Li et al. 2015)
  data_met = generateBetaMatrix(m_met, n)
  if(run == 1 & create.example.plot == TRUE) {
    png(file.path(sub.dir.figures, "Methyl_synth_basic_beta.png"), width = 1000, height = 1000)
    heatmap.2(data_met, dendrogram = "none", trace = "none", Rowv = F,
              Colv = F, labRow = "", main = "Basic Methylation Dataset")
    invisible(dev.off())
  }

  # transform beta value to M values
  if(any(data_met == 1)){
    data_met[which(data_met == 1)] = 1-.Machine$double.eps
  } 
  data_met_M = log2(data_met/(1-data_met))
  
  if(run == 1 & create.example.plot == TRUE) {
    png(file.path(sub.dir.figures, "Methyl_synth_basic.png"), width = 1000, height = 1000)
    heatmap.2(data_met_M, dendrogram = "none", trace = "none", Rowv = F, Colv = F,
              labRow = "", main = "Basic DNA-Methylation Dataset (M-values)")
    invisible(dev.off())
  }  

  # do not add additional random noise to data matrix
  data_met_noise = data_met_M 

  MET_PW_idx = factor(paste0("gr_", sort(rep(c(1:n_groups), m_met/n_groups))))
  saveRDS(MET_PW_idx, file.path(sub.dir.RData, "MET_PWs.RData"))
  
  # Calculate detectable effect in each feature.
  # Effect (diff in means) detectable in t-test of two independent normal distibuted samples 
  # with approximately equal variances. (Sachs p. 385 or Bland p. 336 or Held p. 366) 
  # detectable effect delta
  effect_met = apply(data_met_noise, MARGIN = 1, function(x) {
    power.t.test(delta = NULL, n = (2*n_tu*n_no)/(n_tu+n_no), sig.level = alpha, 
                 power = pow, sd = sd(x), type = "two.sample", alternative = "two.sided")$delta
  })
  
  effect_size_met = effect_met/apply(data_met_noise, MARGIN = 1, sd)
  cat("Relative detectable effect size in each MET feature:", mean(effect_size_met), "\n")
  
  # plot histogram of detectable effect size
  png(file.path(sub.dir.figures, "Effect_size_Methyl.png"))
  hist(effect_met, main = "Effect detectable with n = 200", 
       xlab = "Difference of M-value", breaks = 20)
  invisible(dev.off())
  
  # determine factor for values to add to induce DE
  c_dm = .effect.size
  dx = c_dm/mean(effect_size_met)
  c_dm_abs = dx*median(effect_met)
  cat("Values to add for a median effect of", c_dm, ":", round(c_dm_abs, 3), "\n")

  dx_sd = c_dm_sd/median(effect_size_met)
  c_dm_sd_abs = dx_sd*median(effect_met)
  
  
  if(run == 1 & create.example.plot == TRUE) {
    # plot empirical distribution function of detectable effect size
    # and indicate percentiles of effect to add for DE
    png(file.path(sub.dir.figures, "Methyl_Fn_effect.png"), width = 2000, height = 1000, res = 200)
    plot(ecdf(effect_met), main = "Cumulative distribution of detectable Effect size",
         xlab = "Difference of M-value", xaxt = "n", 
         xlim = c(min(round(min(effect_met),1),  round(c_dm_abs, 1)), 
                  max(round(max(effect_met), 1), round(c_dm_abs, 1))+0.1))
    axis(side = 1, at = c(round(min(effect_met),1), round(c_dm_abs, 3), round(max(effect_met), 1)),
         tick = c(round(min(effect_met),1), round(c_dm_abs, 3), max(round(max(effect_met), 1), max(c_dm))))
    abline(v = c_dm_abs)
    text(x = c_dm_abs, y = 0.8, labels = c_dm, pos = 4)
    invisible(dev.off())
  }

  
  for(k in 1:length(n_groups)){    
    
    if(m_ge %% n_groups[k] == 0 & m_met %% n_groups[k] == 0){

      for(i in 1:length(p_dm_pw)){
        for(j in 1:length(c_dm)){
          
          pw_de_idx = c(sapply(groups_d, function(x) {
            # within perturbed pathways
            sample(which(MET_PW_idx %in% paste0("gr_", x)), size = p_dm_pw[i]/100*m_met/n_groups, 
                   replace = FALSE)
          }))
          
          # randomly spread
          ran_de_idx = sample(1:m_met, size = round(m_met*p_random), replace = FALSE)
          
          diff_idx = unique(c(pw_de_idx, ran_de_idx))
          
          Met_DE = induceDE(data_met_noise, diff_idx, c = c_dm_abs[j], 
                                noise = c_dm_sd_abs, type = "met") #noise = sqrt(0.5))
          
          data_met_de = Met_DE$DE_mtx
          fname = paste("Methyl_synth_", n_groups[k], "gr_", p_dm_pw[i], "p_", pt_dm[j]*100, sep = "")
          saveRDS(data_met_de, file = file.path(sub.dir.RData, paste(fname, "RData", sep = ".")))
          #write.table(data_met_de, file.path(sub.dir.files, paste(fname, "txt", sep = ".")), quote = FALSE)
          if(run == 1 & create.example.plot == TRUE) {
            png(file.path(sub.dir.figures, paste(fname, "png", sep = ".")), width = 1000, height = 1000)
            heatmap.2(data_met_de, dendrogram = "none", trace = "none", Rowv = T, Colv = F,
                      labRow = "", main = paste("\n\nDM Sites in perturbed pathways:", p_dm_pw[i],
                                                "% \n DM Sites in total: ", p_dm[i],
                                                "% \nMean difference between groups in DM Sites:",
                                                "\npercentile: ", pt_dm[j]*100, "/ value: ", c_dm[j]))
            invisible(dev.off())
          }
         
          positives = Met_DE$Diff
          saveRDS(positives, file = file.path(sub.dir.RData, paste0(fname, "_Diff", ".RData")))
          write.table(positives, file.path(sub.dir.files, paste0(fname, "_Diff", ".txt")),
                      row.names = FALSE, quote = FALSE)
          
        }
      }
      
    } else cat("Not all Features assigned to a Group! Quit.\n")
  }

  ttest_result_met = list()
  
  for(k in 1:length(n_groups)){
    for(i in 1: length(p_dm_pw)){
      for(j in 1:length(pt_dm)){
        
        fname = paste("Methyl_synth_", n_groups[k], "gr_", p_dm_pw[i], "p_", pt_dm[j]*100, sep = "")
        data_met_de = readRDS(file.path(sub.dir.RData, paste(fname, "RData", sep = ".")))
        positives = readRDS(file.path(sub.dir.RData, paste0(fname, "_Diff", ".RData")))
        
        ttest_p = apply(data_met_de, 1, function(x) {
          t.test(x[1:n_tu], x[(n_tu+1):n], var.equal = TRUE, paired = FALSE,
                 alternative = "two.sided")$p.value })
        ttest_pFDR = p.adjust(ttest_p, method = "BH")
        positives_t = list(names(which(ttest_pFDR <= alpha)))
        names(positives_t) = fname
        ttest_result_met = c(ttest_result_met, positives_t)
      }
    }
  }
  
  saveRDS(ttest_result_met, file = file.path(sub.dir.RData, "ttest_result_Methyl.RData"))      
      
}

#####------------------------------------------------------------------
# Calculate AUCs and other performance measures for results of t-test 
# on simulated datasets.
#####------------------------------------------------------------------    
source(file.path(.src.dir, "drawROC.R"))

for(run in 1:runs){
  
  # reset sub.dir.paths
  sub.dir.figures = sub.dir.figures.tmp 
  sub.dir.files = sub.dir.files.tmp
  sub.dir.RData = sub.dir.RData.tmp
  
  # set sub.dir.paths to run dir
  sub.dir.figures = file.path(sub.dir.figures, paste0("run", run))
  sub.dir.files = file.path(sub.dir.files, paste0("run", run))
  sub.dir.RData = file.path(sub.dir.RData, paste0("run", run))  
  
  ttest_result_ge = readRDS(file.path(sub.dir.RData, "ttest_result_GeneExp.RData"))
  ttest_result_met = readRDS(file.path(sub.dir.RData, "ttest_result_Methyl.RData"))
  
  ttest_result = c(ttest_result_ge, ttest_result_met)
  
  drawROC(ttest_result, method = "ttest", current.run, stats.run = paste0("run", run))
  cat("Finished drawing ROC curves.\n")
  
  cat("Calculating AUCs, F1-score, Cohen's Kappa and MCC ...\n")
  auc_arr = calculateAUC(ttest_result, method = "ttest", current.run, 
                         stats.run = paste0("run", run))
  
  f1_arr = calculateACC(ttest_result, method = "ttest", current.run,
                       stats.run = paste0("run", run), measure = "f1")
  
  kappa_arr = calculateACC(ttest_result, method = "ttest", current.run,
                        stats.run = paste0("run", run), measure = "kappa")
  
  mcc_arr = calculateACC(ttest_result, method = "ttest", current.run,
                        stats.run = paste0("run", run), measure = "mcc")
  
  
}


cat("Total number of datasets in each stats run: ", 
  length(n_groups)*length(p_de_pw)*length(c_de)*length(.data.to.integrate), 
  "\n")

cat("Finished synthetic data generation.\n")
cat(paste("---------------", Sys.time(), "---------------\n"))


sink()
close(log.con)




