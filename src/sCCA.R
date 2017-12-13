#######################################################################
##### file: sCCA.R                                                #####
##### input: biological or synthetic datasets located in          #####
#####        ./results/.current.run/biologicalData|syntheticData  #####
##### output: vector of penalties and cannonical weights          #####
##### packages: PMA, impute                                       #####
##### author: B. Pucher                                           #####
##### date created: 09/07/2015                                    #####
##### last change:  29/11/2017                                    #####
#######################################################################

rm(list = ls())

#source("http://bioconductor.org/biocLite.R")
#biocLite("impute")

#install.packages("PMA")
library("PMA")

#######################################################################
#####                                                             #####
##### FUNCTIONS                                                   #####
#####                                                             #####
#######################################################################

#####------------------------------------------------------------------
# run sCCA on specified dataset list and save results to .txt and .RData
#####------------------------------------------------------------------
runAndSavesCCA = function(data, penalties, fnames, fname.add = "", 
                          save.txt = TRUE, save.RData = TRUE) {
  
  out = MultiCCA(data, penalty = penalties, ncomponents = 1,
                 niter = .niter, type = "standard", trace = FALSE)
  
  print(out); cat("\n")
  
  sCCA_res = mapply(function(d, o.ws) {
    colnames(d)[which(o.ws != 0)]
  }, d = data, o.ws = out$ws)
  names(sCCA_res) = names(data)
  
  # save selected features to .txt-file
  if(save.txt == TRUE){
    mapply(function(x,n) {
      write.table(x, file.path(sub.dir.files, paste0(n, ".txt")),
                  quote = FALSE, row.names = FALSE, col.names = FALSE)
    }, x = sCCA_res, n = paste0("sCCA_", fnames, fname.add))
    cat("Saved features selected by sCCA to .txt-files. \n")
  }
  
  if(save.RData == TRUE){
    saveRDS(sCCA_res, file.path(sub.dir.RData, paste0("sCCA_result", fname.add, ".RData"))) 
    cat("Saved features selected by sCCA to .RData-files. \n\n")
  }

  return(sCCA_res)
}

#####------------------------------------------------------------------
# Plot distance from desired level of sparsity for the range of 
# penalty terms testd.
#####------------------------------------------------------------------
plotPenalty = function(penalties, dist, fname){
  for(d in 1:nrow(dist)){
    png(file.path(sub.dir.figures, paste0(fname[d], "_penalties.png")), 
        width = 1500, height = 1200, res = 200)
    plot(x = penalties[d,], xlab = "penalty", 
         y = dist[d,], ylim = c(0, max(dist[d,])), ylab = "distance delta n",
         main = paste0(fname[d], ": Distance to desired sparseness"),
         type = "l", lwd = 2, col = "darkblue")
    invisible(dev.off())
  }
}

#####------------------------------------------------------------------
# Create and save array of penalties to test 
#####------------------------------------------------------------------
createPenalties = function(n.penalties, m.features, factor.min, factor.max) {

  cat("Create array with", n.penalties, "sets of penalties to test.. \n")
  penalties <- matrix(NA, nrow = length(.data.to.integrate), ncol = n.penalties)
  for(x in 1:length(.data.to.integrate)){
    # code adappted from function MultiCCA.permute
    penalties[x, ] <- pmax(seq(factor.min, factor.max, len = n.penalties) * 
                             sqrt(m.features[x]), 1.1)
  }
  write.table(round(t(penalties), 2), file.path(sub.dir.files, "penalties_tested.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  saveRDS(penalties, file.path(sub.dir.RData, "penalties_tested.RData"))
  
  return(penalties)
} 

#######################################################################
#####                                                             #####
##### MAIN SECTION                                                #####
#####                                                             #####
#######################################################################
current.run = .current.sCCA
input.data.dir = switch(.type,
                        biolog = .current.biolog,
                        synthet = .current.synthet)
sub.dir.name = file.path(paste0(.type, "sCCA"), current.run)
source(file.path(.src.dir, "setSubDirPath.R"))
source(file.path(.src.dir, "drawROC.R"))
source(file.path(.src.dir, "assessPerformance.R"))

log.con = file(file.path(sub.dir.files, "sCCALog.txt"))
sink(file = log.con, type = "output")
flush(log.con)

cat(paste("---------------", Sys.time(), "---------------"), "\n")
cat("Start calculation of sCCA for a range of tuning parameters ....\n")

# load datasets
data.dir = switch(.type,
                  biolog = file.path("biologicalData", input.data.dir),
                  synthet = file.path("syntheticData", input.data.dir))

datasets.all = list.files(file.path(.run.dir, data.dir, "RData"))
cat("Assuming datasets created in run:", input.data.dir, "\n")

if(.type == "synthet"){
  
  runs = .stats.runs
  stats.seed = .stats.seed
  sub.dir.files.tmp = sub.dir.files
  sub.dir.figures.tmp = sub.dir.figures
  sub.dir.RData.tmp = sub.dir.RData
    
  for(run in 1:runs){
    
    cat("Current stats run:", run, "\n")
    sub.dir.files = file.path(sub.dir.files.tmp, paste0("run", run))
    sub.dir.figures = file.path(sub.dir.figures.tmp, paste0("run", run))
    sub.dir.RData = file.path(sub.dir.RData.tmp, paste0("run", run))
    
    dir.create(sub.dir.files)
    dir.create(sub.dir.figures)
    dir.create(sub.dir.RData)
    
    datasets.all = list.files(file.path(.run.dir, data.dir, "RData", paste0("run", run)))
    datasets = datasets.all[grep("_Diff|train|test|PWs", datasets.all, invert = TRUE)] 
    
    # get n_groups, p_de, p_dm, c_de, c_dm from datasets names
    source(file.path(.src.dir, "helpParameter.R")) 
    
    # create array with penalties to test
    penalties = createPenalties(n.penalties = 250, m.features = c(.m.gene.exp, .m.methylation),
                                factor.min = 0.1, factor.max = 0.6)
    
    sCCA_result = list()  
    num_nz_weights = matrix(NA, nrow = nrow(penalties), ncol = ncol(penalties))
    
    run.seed = stats.seed[run]
    cat("Seed used:", run.seed, "\n")
    set.seed(run.seed)
    
    for(k in 1:length(n_groups)){
      for(i in 1:length(p_de)){
        for(j in 1:length(c_de)){
          
          fname = paste(intersect(c(n_groups[k], "gr_", p_de[i], "p_", c_de[j]), 
                                  c(n_groups[k], "gr_", p_dm[i], "p_", c_dm[j])), collapse = "")
          fname_ge = paste("GeneExp_synth_", n_groups[k], "gr_", p_de[i], "p_", c_de[j], sep = "")
          fname_met = paste("Methyl_synth_",  n_groups[k], "gr_", p_dm[i], "p_", c_dm[j], sep = "")
          data.names = c(fname_ge, fname_met)
          fname = Reduce(intersect, strsplit(data.names, split = "synth"))
          data.paths = as.list(file.path(.run.dir, data.dir, "RData", paste0("run", run),
                                         paste0(data.names, ".RData")))
          data.raw.t = lapply(data.paths, readRDS)
          
          # transpose datasets 
          # MultiCCA expects the same number of rows (samples) in each dataset
          data = lapply(data.raw.t, t); names(data) = data.names
          
          cat("Find penalties for datasets", fname, "\nwith", n_groups[k], "Groups,",
              p_de[i], "% differential profiles and \nan effect of", c_de[j]/100, 
              "between groups.\n\n")
          
          # calculate sparse CCA for penalty pairs and evaluate them
          for(y in 1:ncol(penalties)){
            out = MultiCCA(data, penalty = penalties[,y], 
                           ncomponents = 1, niter = 25, 
                           type = "standard", trace = F)
            
            num_nz_weights[,y] = unlist(lapply(out$ws, function(x) length(x[which(x != 0)])))
          }
          cat("Save number of non-zero elements in canonical weight vectors.\n")
          saveRDS(num_nz_weights, file = file.path(sub.dir.RData, paste0("num_nz_weights", fname, ".RData")))
          
          # determine number of features to be selected by sCCA 
          positives = lapply(file.path(.run.dir, data.dir, "RData", paste0("run", run), 
                                       paste0(c(fname_ge, fname_met), "_Diff.RData")), readRDS)
          positivesL = sapply(positives, length)
          
          to_select = positivesL
          cat("Number of features to be selected by sCCA:", to_select, "\n\n")
          
          dist = abs(num_nz_weights - to_select)
          
          # plot distance to desired sparseness for range of penalties
          plotPenalty(penalties, dist, data.names)
          
          min_idx = apply(dist, MARGIN = 1, which.min)
          bestpenalties = mapply(function(x, y) penalties[x,y], x = c(1:nrow(penalties)), y = min_idx)
          
          current_res = runAndSavesCCA(data, penalties = bestpenalties, fnames = data.names, 
                                       save.RData = FALSE)
          
          sCCA_result = c(sCCA_result, current_res)
          
          # calculate statistics for penalty terms
          perm.out = MultiCCA.permute(data, penalties = as.matrix(bestpenalties, nrow = 2),
                                      type = "standard", nperms = 1000, trace = FALSE)
          
          print(perm.out); cat("\n")
          
        }
      }
    }
    
    cat("Finished calculation of MultiCCA.\n")
    cat(paste("---------------", Sys.time(), "---------------"), "\n\n")
    
    save(sCCA_result, file = file.path(sub.dir.RData, "sCCA_result.RData"))

    drawROC(sCCA_result, method = "sCCA", input.data.dir, stats.run = paste0("run", run))
    cat("Finished drawing ROC curves.\n")
 
    cat("Calculating AUCs, ACC, F1-score, Cohen's Kappa and MCC ...\n")   
    auc_arr = calculateAUC(sCCA_result, method = "sCCA", input.data.dir, 
                           stats.run = paste0("run", run))
  
    acc_arr = calculateACC(sCCA_result, method = "sCCA", input.data.dir, 
                           stats.run = paste0("run", run), measure = "acc")
    
    f1_arr = calculateACC(sCCA_result, method = "sCCA", input.data.dir,
                         stats.run = paste0("run", run), measure = "f1")
    
    kappa_arr = calculateACC(sCCA_result, method = "sCCA", input.data.dir,
                          stats.run = paste0("run", run), measure = "kappa")
    
    mcc_arr = calculateACC(sCCA_result, method = "sCCA", input.data.dir,
                             stats.run = paste0("run", run), measure = "mcc")
    
  }
  
  # store results of stats runs in one list
  sCCA_result_list = vector("list", runs) 
  for(run in 1:runs){
    
    sub.dir.RData = file.path(sub.dir.RData.tmp, paste0("run", run))
    load(file.path(sub.dir.RData, "sCCA_result.RData"))
    sCCA_result_list[[run]] = sCCA_result
  }
  
  saveRDS(sCCA_result_list, file = file.path(sub.dir.RData.tmp, "sCCA_result_list.RData"))

} else {
  
  datasets = datasets.all[grep("train|test", datasets.all, invert = TRUE)]
  data.names = gsub(".RData", "", datasets)
  data.paths = as.list(file.path(.run.dir, data.dir, "RData", datasets))
  data.raw.t = lapply(data.paths, readRDS)

  # transpose datasets 
  # MultiCCA expects the same number of rows (samples) in each dataset
  data = lapply(data.raw.t, t)
  names(data) = data.names
  
  # create array of penalties to test 
  penalties = createPenalties(n.penalties = 100, m.features = sapply(data, ncol),
                              factor.min = 0.01, factor.max = 0.1)

  # calculate sparse CCA for penalty pairs 
  num_nz_weights = matrix(NA, nrow = nrow(penalties), ncol = ncol(penalties))
  for(y in 1:ncol(penalties)){
    out = MultiCCA(data, penalty = penalties[,y], 
                   ncomponents = 1, niter = 25, 
                   type = "standard", trace = FALSE)
    
    num_nz_weights[,y] = unlist(lapply(out$ws, function(x) length(x[which(x != 0)]))) 
  }
  
  cat("Save number of non-zero elements in canonical weight vectors.\n")
  saveRDS(num_nz_weights, file = file.path(sub.dir.RData, "num_nz_weights.RData"))

  cat("Finished calculation of sCCA for range of tuning parameters.\n")
  cat(paste("-----", Sys.time(), "-----"), "\n\n")
  cat("Extract penalties for desired level of sparseness.\n")

  # determine number of features to be selected by sCCA 
  to_select = round(unlist(lapply(lapply(data, dim), "[", 2)) * (.pr/100))
  cat("Number of features to be selected by sCCA:", to_select, "\n\n")
  
  dist = abs(num_nz_weights - to_select)
  
  # plot distance to desired sparseness for range of penalties
  plotPenalty(penalties, dist, data.names)

  min_idx = apply(dist, MARGIN = 1, which.min)
  bestpenalties = mapply(function(x, y) penalties[x,y], x = c(1:nrow(penalties)), y = min_idx)

  cat("Starting cross validation ...\n")
  cat(paste("-----", Sys.time(), "-----"), "\n\n")
  
  # do cross validation
  if(.do.CV == TRUE) {
    
    source(file.path(.src.dir, "crossValidation.R"))
    
    # run sCCA on training sets
    for(s in 1:.subsets.val){
      data.path.train = file.path(.run.dir, data.dir, "RData", paste0("data_train", s, ".RData"))
      data.train.t = readRDS(data.path.train)
      
      # transpose datasets
      data.train = lapply(data.train.t, t)
      
      cat("Run sCCA on training set", s, ": \n")
      sCCA_result_train = runAndSavesCCA(data.train, bestpenalties, fnames = names(data.train), 
                                         fname.add = paste0("_train", s))
      
      # calculate statistics for penalty terms
      set.seed(.biolog.seed.cv[s])
      cat("Calculate statistics for penalty terms:\nSeed used:", .biolog.seed.cv[s], "\n")
      perm.out = MultiCCA.permute(data.train, penalties = as.matrix(bestpenalties, nrow = 2),
                                  type = "standard", nperms = 1000, trace = FALSE)
      print(perm.out)

    }
    
    # evaluate performance
    evaluatePerformance(method = "sCCA", subsets = .subsets.val)

  }
  
  cat(paste("-----", Sys.time(), "-----"), "\n\n")
  
  # run sCCA on biological data
  sCCA_result = runAndSavesCCA(data, bestpenalties, fnames = data.names)
  
  # calculate statistics for penalty terms
  set.seed(.biolog.seed)
  perm.out = MultiCCA.permute(data, penalties = as.matrix(bestpenalties, nrow = 2),
                              type = "standard", nperms = 1000, trace = FALSE)
  print(perm.out)
  
  cat(paste("-----", Sys.time(), "-----"), "\n\n")
  cat("\nFinished calculation of MultiCCA.\n")
  
}

sink()
close(log.con)

