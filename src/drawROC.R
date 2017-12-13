#######################################################################
##### file: drawROC.R                                             #####
##### description: calculation of TPR and FPR, AUC and other      #####
#####   performance measures; construction of ROC curves          #####
##### input: list containing the genes identified by method       #####
##### output: TPR, FPR, ROC curve plots                           #####
##### author: B. Pucher                                           #####
##### date created: 16/03/2017                                    #####
##### last change: 09/11/2017                                     #####
#######################################################################

#####------------------------------------------------------------------
# Calculate true positive rate (TPR) and false positive rate (FPR)
# for a given feature list resulting from method and draw ROC curves.
#####------------------------------------------------------------------
drawROC = function(feature_list, method, synthet.run, stats.run = ""){
  
  # For each dataset: grep the features detected by a method from
  # feature_list, calculate confusion matrix and TPR and FPR.
  # TPR and FPR for all datasets are saved to an array.
  # Based on array of TPRs and FPRs draw ROC curves.
  
  datasets = paste0(names(feature_list), ".RData")

  n_groups = sort(as.numeric(unlist(unique(strsplit(unlist(lapply(strsplit(datasets, split = "_"), "[", 3)), split = "gr")))))
  p_de = sort(as.numeric(unlist(unique(strsplit(unlist(lapply(strsplit(datasets[grep("GeneExp", datasets)], split = "_"), "[", 4)), split = "p")))))
  p_dm = sort(as.numeric(unlist(unique(strsplit(unlist(lapply(strsplit(datasets[grep("Methyl", datasets)], split = "_"), "[", 4)), split = "p")))))
  c_de = sort(as.numeric(unlist(unique(strsplit(unlist(lapply(strsplit(datasets[grep("GeneExp", datasets)], split = "_"), "[", 5)), split = ".RData")))))
  c_dm = sort(as.numeric(unlist(unique(strsplit(unlist(lapply(strsplit(datasets[grep("Methyl", datasets)], split = "_"), "[", 5)), split = ".RData")))))
  
  TPR_arr = array(dim = c(length(p_de), length(c_de), length(n_groups)))
  dimnames(TPR_arr) = list(paste0(p_de, "%"), paste0("+", c_de), paste(n_groups, "Gr"))
  FPR_arr = array(dim = c(length(p_dm), length(c_dm), length(n_groups)))
  dimnames(FPR_arr) = list(paste0(p_de, "%"), paste0("+", c_de), paste(n_groups, "Gr"))
  
  for(k in 1:length(n_groups)){
    for(i in 1:length(p_de)){
      for(j in 1:length(c_de)){
        
        fname_ge = paste("GeneExp_synth_", n_groups[k], "gr_", p_de[i], "p_", c_de[j], sep = "")
        fname_met = paste("Methyl_synth_",  n_groups[k], "gr_", p_dm[i], "p_", c_dm[j], sep = "")
        data.names = c(fname_ge, fname_met)
        
        # POSITIVES: differentially expressed/methylated features
        pos_ge = readRDS(file.path(.run.dir, "syntheticData", synthet.run, "RData", 
                                   stats.run, paste0(fname_ge, "_Diff.RData")))
        pos_met = readRDS(file.path(.run.dir, "syntheticData", synthet.run, "RData", 
                                    stats.run, paste0(fname_met, "_Diff.RData")))
        
        data_tmp_ge = readRDS(file.path(.run.dir, "syntheticData", synthet.run, "RData", 
                                        stats.run, paste0(fname_ge, ".RData")))
        data_tmp_met = readRDS(file.path(.run.dir, "syntheticData", synthet.run, "RData", 
                                         stats.run, paste0(fname_met, ".RData")))
        
        # calculate confusion matrices for current datasets
        cat("Confusion matrices for", data.names, "\n")
        
        TP_ge = intersect(pos_ge, feature_list[[fname_ge]]); cat("TP genes:", length(TP_ge), "\n")
        FP_ge = setdiff(feature_list[[fname_ge]], pos_ge); cat("FP genes:", length(FP_ge[which(FP_ge != "")]), "\n")
        FN_ge = setdiff(pos_ge, feature_list[[fname_ge]]); cat("FN genes:", length(FN_ge[which(FN_ge != "")]), "\n")
        TN_ge = setdiff(rownames(data_tmp_ge), c(TP_ge, FP_ge, FN_ge)); cat("TN genes:", length(TN_ge[which(TN_ge != "")]), "\n")
        
        TP_met = intersect(pos_met, feature_list[[fname_met]]); cat("TP methyl:", length(TP_met), "\n")
        FP_met = setdiff(feature_list[[fname_met]], pos_met); cat("FP methyl:", length(FP_met[which(FP_met != "")]), "\n")
        FN_met = setdiff(pos_met, feature_list[[fname_met]]); cat("FN methyl:", length(FN_met[which(FN_met != "")]), "\n")
        TN_met = setdiff(rownames(data_tmp_met), c(TP_met, FP_met, FN_met)); cat("TN methyl:", length(TN_met[which(TN_met != "")]), "\n")
        
        TPR = (length(TP_ge)+length(TP_met))/(length(pos_ge)+length(pos_met))
        FPR = (length(FP_ge[which(FP_ge != "")])+length(FP_met[which(FP_met != "")]))/
          (length(TN_ge[which(TN_ge != "")])+length(TN_met[which(TN_met != "")])+
             length(FP_ge[which(FP_ge != "")])+length(FP_met[which(FP_met != "")]))
        cat("TPR:", TPR, "\n")
        cat("FPR:", FPR, "\n")
        cat("\n\n")
        
        TPR_arr[i,j,k] = TPR; 
        FPR_arr[i,j,k] = FPR; 
      }
    }
  }
  
  saveRDS(TPR_arr, file = file.path(sub.dir.RData, paste0(method, "_TPR_arr.RData")))
  saveRDS(FPR_arr, file = file.path(sub.dir.RData, paste0(method, "_FPR_arr.RData")))
  
  cat("Start plotting ROC curves...\n")
  
  # Plot ROC curve for % of differential profiles
  for(k in 1:length(n_groups)){
      p_de_total = 1/n_groups[k]*(n_groups[k]*.p.groups.d)*.percentage.pw + .p.random*100
      png(file.path(sub.dir.figures, paste0(method, "_ROC_", n_groups[k], "gr_diff_value.png")))
    plot(1, type="n", axes=F, xlab="", ylab="") 
    par(new = FALSE)
    for(j in 1:length(c_de)){
      plot(FPR_arr[,j,k], TPR_arr[,j,k], xlim = c(-0.1,1.1), ylim = c(-0.1,1.1), 
           main = paste("ROC of", method, ":", n_groups[k], "Groups"), xlab = "FPR", ylab = "TPR", 
           col = rainbow(length(c_de))[j], lwd = 3, type = "b")
      text(FPR_arr[,j,k], TPR_arr[,j,k], labels = paste(p_de_total, "%"), pos = 2)
      par(new = TRUE)
    }
    legend('bottomright', paste("effect size:", c_de, "%"), col = rainbow(length(c_de)),
           lty = 1, cex = 1)
    dev.off()
  }
  cat("ROC comparing percentage of differential profiles plotted. \n")
  
  # Plot ROC curve for differences between conditions
  for(k in 1:length(n_groups)){
    p_de_total = 1/n_groups[k]*(n_groups[k]*.p.groups.d)*.percentage.pw + .p.random*100
    png(file.path(sub.dir.figures, paste0(method, "_ROC_", n_groups[k], "gr_diff_percentage.png"))) 
    par(new = FALSE)
    for(i in 1:length(p_de)){
      plot(FPR_arr[i,,k], TPR_arr[i,,k], xlim = c(-0.1,1.1), ylim = c(-0.1,1.1), 
           main = paste("ROC of", method, ":", n_groups[k], "Groups"),
           xlab = "FPR", ylab = "TPR", col = rainbow(length(p_de_total))[i], lwd = 3, type = "b")
      text(FPR_arr[i,,k], TPR_arr[i,,k], labels = paste0(c_de, "%"), pos = 4)
      par(new = TRUE)
    }
    legend('bottomright', paste(p_de_total, "% differential profiles"), 
           col = rainbow(length(p_de_total)), lty = 1, cex = 1)
    dev.off()
  }
  cat("ROC comparing mean difference of differential profiles between conditions plotted. \n")
  
}

#####------------------------------------------------------------------
# Calculate the area under the curve with package pROC
#####------------------------------------------------------------------
library(pROC)
calculateAUC = function(feature_list, method, synthet.run, stats.run){
  
  datasets = paste0(names(feature_list), ".RData")
  
  n_groups = sort(as.numeric(unlist(unique(strsplit(unlist(lapply(strsplit(datasets, split = "_"), "[", 3)), split = "gr")))))
  p_de = sort(as.numeric(unlist(unique(strsplit(unlist(lapply(strsplit(datasets[grep("GeneExp", datasets)], split = "_"), "[", 4)), split = "p")))))
  p_dm = sort(as.numeric(unlist(unique(strsplit(unlist(lapply(strsplit(datasets[grep("Methyl", datasets)], split = "_"), "[", 4)), split = "p")))))
  c_de = sort(as.numeric(unlist(unique(strsplit(unlist(lapply(strsplit(datasets[grep("GeneExp", datasets)], split = "_"), "[", 5)), split = ".RData")))))
  c_dm = sort(as.numeric(unlist(unique(strsplit(unlist(lapply(strsplit(datasets[grep("Methyl", datasets)], split = "_"), "[", 5)), split = ".RData")))))
  
  auc_arr = array(dim = c(length(p_de), length(c_de), length(n_groups)))
  dimnames(auc_arr) = list(paste0(p_de, "%"), paste0("+", c_de), paste(n_groups, "Gr"))
  
  for(k in 1:length(n_groups)){
    for(i in 1:length(p_de)){
      for(j in 1:length(c_de)){
        
        fname_ge = paste("GeneExp_synth_", n_groups[k], "gr_", p_de[i], "p_", c_de[j], sep = "")
        fname_met = paste("Methyl_synth_",  n_groups[k], "gr_", p_dm[i], "p_", c_dm[j], sep = "")
        data.names = c(fname_ge, fname_met)
        
        data_tmp_ge = readRDS(file.path(.run.dir, "syntheticData", synthet.run, "RData", 
                                        stats.run, paste0(fname_ge, ".RData")))
        data_tmp_met = readRDS(file.path(.run.dir, "syntheticData", synthet.run, "RData", 
                                         stats.run, paste0(fname_met, ".RData")))
        
        # RESPONSE: differentially expressed/methylated features
        pos_ge = readRDS(file.path(.run.dir, "syntheticData", synthet.run, "RData", 
                                   stats.run, paste0(fname_ge, "_Diff.RData")))
        pos_met = readRDS(file.path(.run.dir, "syntheticData", synthet.run, "RData", 
                                    stats.run, paste0(fname_met, "_Diff.RData")))
        resp_ge = 1* (rownames(data_tmp_ge) %in% pos_ge)
        resp_met = 1*(rownames(data_tmp_met) %in% pos_met)
        resp = c(resp_ge, resp_met)
          
        # PREDICTOR: features selected by method
        pred_ge = 1*(rownames(data_tmp_ge) %in% feature_list[[fname_ge]])
        pred_met = 1*(rownames(data_tmp_met) %in% feature_list[[fname_met]])
        pred = c(pred_ge, pred_met)
        
        auc_arr[i,j,k] = auc(resp, pred)
      }
    }
  }

  saveRDS(auc_arr, file = file.path(sub.dir.RData, paste0(method, "_auc_arr.RData")))
  cat("Saved AUC of ROC curves to RData folder.\n")
  
  return(auc_arr)
}

#####------------------------------------------------------------------
# Calculate the accuracy of methods in identifying DE features 
#####------------------------------------------------------------------
calculateACC = function(feature_list, method, synthet.run, stats.run,
                        measure = "acc"){
  
  # parameter help: 
  # measure = c("acc", "f1", "kappa", "mcc")
  
  datasets = paste0(names(feature_list), ".RData")
  
  n_groups = sort(as.numeric(unlist(unique(strsplit(unlist(lapply(strsplit(datasets, split = "_"), "[", 3)), split = "gr")))))
  p_de = sort(as.numeric(unlist(unique(strsplit(unlist(lapply(strsplit(datasets[grep("GeneExp", datasets)], split = "_"), "[", 4)), split = "p")))))
  p_dm = sort(as.numeric(unlist(unique(strsplit(unlist(lapply(strsplit(datasets[grep("Methyl", datasets)], split = "_"), "[", 4)), split = "p")))))
  c_de = sort(as.numeric(unlist(unique(strsplit(unlist(lapply(strsplit(datasets[grep("GeneExp", datasets)], split = "_"), "[", 5)), split = ".RData")))))
  c_dm = sort(as.numeric(unlist(unique(strsplit(unlist(lapply(strsplit(datasets[grep("Methyl", datasets)], split = "_"), "[", 5)), split = ".RData")))))
  
  acc_arr = array(dim = c(length(p_de), length(c_de), length(n_groups)))
  dimnames(acc_arr) = list(paste0(p_de, "%"), paste0("+", c_de), paste(n_groups, "Gr"))
  
  for(k in 1:length(n_groups)){
    for(i in 1:length(p_de)){
      for(j in 1:length(c_de)){
        
        fname_ge = paste("GeneExp_synth_", n_groups[k], "gr_", p_de[i], "p_", c_de[j], sep = "")
        fname_met = paste("Methyl_synth_",  n_groups[k], "gr_", p_dm[i], "p_", c_dm[j], sep = "")
        data.names = c(fname_ge, fname_met)
        
        data_tmp_ge = readRDS(file.path(.run.dir, "syntheticData", synthet.run, "RData", 
                                        stats.run, paste0(fname_ge, ".RData")))
        data_tmp_met = readRDS(file.path(.run.dir, "syntheticData", synthet.run, "RData", 
                                         stats.run, paste0(fname_met, ".RData")))
        
        pos_ge = readRDS(file.path(.run.dir, "syntheticData", synthet.run, "RData", 
                                   stats.run, paste0(fname_ge, "_Diff.RData")))
        pos_met = readRDS(file.path(.run.dir, "syntheticData", synthet.run, "RData", 
                                    stats.run, paste0(fname_met, "_Diff.RData")))
        
        TP_ge = length(intersect(pos_ge, feature_list[[fname_ge]]))
        FP_ge = length(setdiff(feature_list[[fname_ge]], pos_ge))
        FN_ge = length(setdiff(pos_ge, feature_list[[fname_ge]]))
        TN_ge = length(rownames(data_tmp_ge)) - (TP_ge + FP_ge + FN_ge)
        
        TP_met = length(intersect(pos_met, feature_list[[fname_met]]))
        FP_met = length(setdiff(feature_list[[fname_met]], pos_met))
        FN_met = length(setdiff(pos_met, feature_list[[fname_met]]))
        TN_met = length(rownames(data_tmp_met)) - (TP_met+FP_met+FN_met) 
        
        TP = as.numeric(TP_ge + TP_met) # force to be numeric not integer
        FP = as.numeric(FP_ge + FP_met)
        TN = as.numeric(TN_ge + TN_met)
        FN = as.numeric(FN_ge + FN_met)
        
        base_table = matrix(c(TP, FP, FN, TN), nrow = 2, byrow = TRUE)
        
        po = sum(diag(base_table))/sum(base_table)
        pe = sum(colSums(t(base_table))*colSums(base_table))/sum(base_table)^2
        
        acc_arr[i,j,k] = switch(measure,
                                # accuracy
                                acc = (TP+TN)/(TP+FN+TN+FP),
                                
                                # F1 score
                                f1 = (2*TP)/(2*TP+FP+FN),
                                
                                # Kappa statistic
                                kappa = (po-pe)/(1-pe),
                                
                                # Matthew's correlation coeffiecient
                                mcc = (TP*TN - FP*FN)/
                                  sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
                                )
                                
                                
                            
      }
    }
  }
  
  saveRDS(acc_arr, file = file.path(sub.dir.RData, paste(method, measure, "arr.RData", sep = "_")))
  cat("Saved", toupper(measure), "of methods on simulated datasets to RData folder.\n")
  
  return(acc_arr)
}

#####------------------------------------------------------------------
# Calculate and save base table with TP, FP, TN and FN.  
#####------------------------------------------------------------------
calculateBaseTab = function(feature_list, method, synthet.run, stats.run){
  
  datasets = paste0(names(feature_list), ".RData")
  
  n_groups = sort(as.numeric(unlist(unique(strsplit(unlist(lapply(strsplit(datasets, split = "_"), "[", 3)), split = "gr")))))
  p_de = sort(as.numeric(unlist(unique(strsplit(unlist(lapply(strsplit(datasets[grep("GeneExp", datasets)], split = "_"), "[", 4)), split = "p")))))
  p_dm = sort(as.numeric(unlist(unique(strsplit(unlist(lapply(strsplit(datasets[grep("Methyl", datasets)], split = "_"), "[", 4)), split = "p")))))
  c_de = sort(as.numeric(unlist(unique(strsplit(unlist(lapply(strsplit(datasets[grep("GeneExp", datasets)], split = "_"), "[", 5)), split = ".RData")))))
  c_dm = sort(as.numeric(unlist(unique(strsplit(unlist(lapply(strsplit(datasets[grep("Methyl", datasets)], split = "_"), "[", 5)), split = ".RData")))))
  
  base_tabs = list()
  
  for(k in 1:length(n_groups)){
    for(i in 1:length(p_de)){
      for(j in 1:length(c_de)){
        
        fname_ge = paste("GeneExp_synth_", n_groups[k], "gr_", p_de[i], "p_", c_de[j], sep = "")
        fname_met = paste("Methyl_synth_",  n_groups[k], "gr_", p_dm[i], "p_", c_dm[j], sep = "")
        data.names = c(fname_ge, fname_met)
        
        data_tmp_ge = readRDS(file.path(.run.dir, "syntheticData", synthet.run, "RData", 
                                        stats.run, paste0(fname_ge, ".RData")))
        data_tmp_met = readRDS(file.path(.run.dir, "syntheticData", synthet.run, "RData", 
                                         stats.run, paste0(fname_met, ".RData")))
        
        pos_ge = readRDS(file.path(.run.dir, "syntheticData", synthet.run, "RData", 
                                   stats.run, paste0(fname_ge, "_Diff.RData")))
        pos_met = readRDS(file.path(.run.dir, "syntheticData", synthet.run, "RData", 
                                    stats.run, paste0(fname_met, "_Diff.RData")))
        
        TP_ge = length(intersect(pos_ge, feature_list[[fname_ge]]))
        FP_ge = length(setdiff(feature_list[[fname_ge]], pos_ge))
        FN_ge = length(setdiff(pos_ge, feature_list[[fname_ge]]))
        TN_ge = length(rownames(data_tmp_ge)) - (TP_ge + FP_ge + FN_ge)
        
        base_table_ge = matrix(c(TP_ge, FP_ge, FN_ge, TN_ge),
                               nrow = 2, byrow = TRUE, 
                               dimnames = list(c("PR_pos", "PR_neg"), 
                                               c("RE_pos", "RE_neg")))
        
        TP_met = length(intersect(pos_met, feature_list[[fname_met]]))
        FP_met = length(setdiff(feature_list[[fname_met]], pos_met))
        FN_met = length(setdiff(pos_met, feature_list[[fname_met]]))
        TN_met = length(rownames(data_tmp_met)) - (TP_met+FP_met+FN_met) 
        
        base_table_met = matrix(c(TP_met, FP_met, FN_met, TN_met),
                               nrow = 2, byrow = TRUE,
                               dimnames = list(c("PR_pos", "PR_neg"), 
                                               c("RE_pos", "RE_neg")))
        
        base_table = list(base_table_ge, base_table_met)
        names(base_table) = data.names
        base_tabs = c(base_tabs, base_table)
        
      }
    }
  }
  
  saveRDS(base_tabs, file = file.path(sub.dir.RData, paste0(method, "_base_tables.RData")))
  cat("Saved Base Tables of methods on simulated datasets to RData folder.\n")
  
  return(base_tabs)
}

#####------------------------------------------------------------------
# Calculate AUC for pathway over-representation analysis.
#####------------------------------------------------------------------
calculatePwAUC = function(feature_list, PW_idx, method, synthet.run, 
                          stats.run, dest.dir, measure = "mcc"){
  
  # parameter help: measure = c("auc", "mcc")
 
  datasets = paste0(names(feature_list), ".RData")
  
  # get n_groups, p_de, p_dm, c_de, c_dm from datasets names
  n_groups = sort(as.numeric(unlist(unique(strsplit(unlist(lapply(strsplit(datasets, split = "_"), "[", 3)), split = "gr")))))
  p_de = sort(as.numeric(unlist(unique(strsplit(unlist(lapply(strsplit(datasets[grep("GeneExp", datasets)], split = "_"), "[", 4)), split = "p")))))
  p_dm = sort(as.numeric(unlist(unique(strsplit(unlist(lapply(strsplit(datasets[grep("Methyl", datasets)], split = "_"), "[", 4)), split = "p")))))
  c_de = sort(as.numeric(unlist(unique(strsplit(unlist(lapply(strsplit(datasets[grep("GeneExp", datasets)], split = "_"), "[", 5)), split = ".RData")))))
  c_dm = sort(as.numeric(unlist(unique(strsplit(unlist(lapply(strsplit(datasets[grep("Methyl", datasets)], split = "_"), "[", 5)), split = ".RData")))))
  
  auc_arr = array(dim = c(length(p_de), length(c_de), length(n_groups)))
  dimnames(auc_arr) = list(paste0(p_de, "%"), paste0("+", c_de), paste(n_groups, "Gr"))
  
  # number of perturbed pathways
  n_groups_d = round(.n.groups*.p.groups.d)
  
  for(k in 1:length(n_groups)){
    for(i in 1:length(p_de)){
      for(j in 1:length(c_de)){
        
        fname_ge = paste("GeneExp_synth_", n_groups[k], "gr_", p_de[i], "p_", c_de[j], sep = "")
        fname_met = paste("Methyl_synth_",  n_groups[k], "gr_", p_dm[i], "p_", c_dm[j], sep = "")
        data.names = c(fname_ge, fname_met)
        add = Reduce(intersect, strsplit(data.names, split = "synth_"))
        data.paths = as.list(file.path(.run.dir, data.dir, "RData", paste0("run", stats.run),
                                       paste0(data.names, ".RData")))
        data.raw = lapply(data.paths, readRDS)
        names(data.raw) = data.names
        
        # determine which features are members of which pathway
        PW_fts = lapply(unique(unlist(PW_idx)), function(pw) {
          fts = mapply(function(x, y) {
            rownames(x)[which(y == pw)]
          }, x = data.raw, y = PW_idx) })
        names(PW_fts) = as.character(unique(unlist(PW_idx)))
        
        ft.universe = sapply(data.raw, rownames)
        current.res = list(GeneExp = feature_list[[fname_ge]],
                           Methyl = feature_list[[fname_met]])
        
        # do over-representation analysis
        cat(data.names, "\n")
        PW_p_temp = lapply(lapply(PW_fts, unlist), testOverrep, 
                        features = unlist(current.res),
                        universe = unlist(ft.universe))
        names(PW_p_temp) <- names(PW_fts)
        PW_p = unlist(PW_p_temp, use.names = TRUE)
          
        cat(method, ": Pathways with adjusted p-values < 0.05 in right-tailed Fisher's test:\n",
              "(", length(PW_p[which(p.adjust(PW_p, method = "BH") < 0.05)]), "):",
              PW_p[which(p.adjust(PW_p, method = "BH") < 0.05)], "\n")
        
        PW_sign = names(PW_p)[which(p.adjust(PW_p, method = "BH") < 0.05)]
        
        # pathways perturbed in current stats run
        set.seed(.stats.seed[stats.run])
        groups_d = sample(n_groups[k], n_groups_d[k], replace = FALSE)
        PW_d = paste0("gr_", groups_d)
          
        resp = 1*names(PW_fts) %in% PW_d
        pred = 1*names(PW_fts) %in% PW_sign
        
        TP = as.numeric(length(intersect(PW_sign, PW_d)))
        FP = as.numeric(length(setdiff(PW_sign, PW_d)))
        FN = as.numeric(length(setdiff(PW_d, PW_sign)))
        TN = as.numeric(length(names(PW_fts)) - (TP+FP+FN))
        
        mcc_denom = ifelse((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) == 0, yes = sqrt(1), 
                           no = sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)))
          
        auc_arr[i,j,k] = switch(measure,
                                
                                # area under the ROC curve
                                auc = auc(resp, pred),
                                
                                # Matthews correlation coefficient
                                mcc = (TP*TN - FP*FN)/mcc_denom
                                )
        
        
      }
    }
  }
  
  saveRDS(auc_arr, file = file.path(dest.dir, paste0(method, "_PW_", measure, "_arr.RData")))
  cat("Saved Pathway", toupper(measure), "to RData folder.\n")
}

#####------------------------------------------------------------------
# Over-representation analysis of pathways in synthetic data.  
#####------------------------------------------------------------------
testOverrep = function(features, universe, pathway){
  
  a = intersect(features, pathway)
  b = setdiff(pathway, features)
  c = intersect(setdiff(universe, pathway), features)
  d = intersect(setdiff(universe, pathway), setdiff(universe, features))
  
  ct =  matrix(c(length(a), length(b), length(c), length(d)),
               nrow = 2, byrow = FALSE,
               dimnames = list("Pathway" = c("Yes", "No"), 
                               "Method" = c("Yes", "No")))
  fp = fisher.test(ct, alternative = "greater")$p.value

  return(fp)
}
      




