#######################################################################
##### file: 5c_postprocOfNMF.R                                    #####
##### input: output of NMF: best factor matrices W and H          #####
##### output: flat list of candidate genes                        #####
##### packages: --                                                #####
##### author: B. Pucher                                           #####
##### date created: 13/04/2015                                    #####
##### last change: 29/11/2017                                     #####
#######################################################################
rm(list = ls())

#######################################################################
#####                                                             #####
##### FUNCTIONS                                                   #####
#####                                                             #####
#######################################################################

#####------------------------------------------------------------------
# Get desired nb of features to be identified by NMF in biologic data
#####------------------------------------------------------------------
getDesiredBiolog = function(bestWHlist){
  
  data.names = bestWHlist$datasets

  # Number of simulated POSITIVES:
  desired = mapply(function(h, p) {
    round((dim(h)[2] * (p/100))/2)
  }, h = bestWHlist[grep("H", names(bestWHlist))], p = .pr, SIMPLIFY = FALSE)
    
  names(desired) = data.names
  
  return(desired)
}

#####------------------------------------------------------------------
# Get desired nb of features to be identified by NMF in synthetic data
#####------------------------------------------------------------------
getDesiredSynthet = function(bestWHlist, synthet.run = "", 
                             stats.run = ""){
  
  data.names = bestWHlist$datasets
  fname_ge = data.names[grep("GeneExp", data.names)]
  fname_met = data.names[grep("Methyl", data.names)]
  
  # POSITIVES: differentially expressed/methylated features
  pos_ge = readRDS(file.path(.run.dir, "syntheticData", synthet.run, 
                             "RData", stats.run, paste0(fname_ge, "_Diff.RData")))
  pos_met = readRDS(file.path(.run.dir, "syntheticData", synthet.run, 
                              "RData", stats.run, paste0(fname_met, "_Diff.RData")))
  
  desired = list(length(pos_ge), length(pos_met))

  return(desired)
}

#####------------------------------------------------------------------
# Transformation of feature weights to z-values
#####------------------------------------------------------------------
zTransform = function(mtx, mode = "mean") {
  
  switch(mode,
         mean = {
           # calculate row-wise mean and standard deviation 
           M = apply(mtx, MARGIN = 1, FUN = mean)
           V = apply(mtx, MARGIN = 1, FUN = sd)
         },
         median = { 
           # calculate row-wise median and median standard deviation 
           M = apply(mtx, MARGIN = 1, FUN = mode)
           V = apply(mtx, MARGIN = 1, FUN = mad, constant = 1)
         })
  
  # transformation of weights to z-values 
  mtx_z = mapply(function(r, m, v) {
    (mtx[r,] - m)/v
  }, r = 1:nrow(mtx), M, V)
  
  return(t(mtx_z))
}

#####------------------------------------------------------------------
# Extract NMF features from WHlist and save results to .txt and .RData
#####------------------------------------------------------------------
extractAndSave = function(WHlist, desired, fnames, fname.add = "",
                          fname.add2 = ""){
  
  W = WHlist$W
  Hi = WHlist[grep("H", names(WHlist))]
  Hi_z = lapply(Hi, zTransform, mode = "median")
  
  good.k = dim(W)[2]
  
  NMF_res = mapply(function(h, d) {
    # concatenate weights of features in a data type for all md-modules 
    # preserving names
    h.vc.list = lapply(1:nrow(h), function(x) { structure(h[x,], names = dimnames(h)[[2]]) })
    h.vc = do.call(c, h.vc.list)
    # order weights 
    h.sort = sort(h.vc, decreasing = TRUE)
    names(h.sort) = sub(".", "", names(h.sort), fixed = TRUE) 
    # remove duplicated feature names from ordered weight vector to get 
    # top unique features
    h.u = h.sort[-which(duplicated(names(h.sort)))]
    names(h.u)[1:d]
  }, h = Hi_z, d = desired)
  names(NMF_res) = WHlist$datasets   
  
  invisible(mapply(function(x, n) {
    write.table(x, file.path(sub.dir.files, paste0(n, ".txt")),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
  }, x = NMF_res, n = paste0("NMF_", fnames, fname.add2)))
  cat("Saved features selected by NMF to .txt-files. \n")
  
  saveRDS(NMF_res, file = file.path(sub.dir.RData, paste0("NMF_result", fname.add, 
                                                fname.add2, ".RData")))
  cat("Saved features selected by NMF to .RData-files. \n\n")
  
  return(NMF_res)
}

#######################################################################
#####                                                             #####
##### MAIN SECTION                                                #####
#####                                                             #####
#######################################################################
current.run = .current.NMF
input.data.dir = switch(.type,
                        biolog = .current.biolog,
                        synthet = .current.synthet)
sub.dir.name = file.path(paste0(.type, "NMF"), current.run)
source(file.path(.src.dir, "setSubDirPath.R"))
source(file.path(.src.dir, "drawROC.R"))
source(file.path(.src.dir, "assessPerformance.R"))

log.con = file(
  file.path(sub.dir.files, paste0("postprocOfNMFLog_", Sys.Date(), ".txt"))
)
sink(file = log.con, type = "output")
flush(log.con)

cat("Start extracting features detected by NMF ...\n")
cat("Assuming datasets from run:", input.data.dir, "\n")
cat("Grep number of features in matrix H of factorization according to mode.\n")

# load datasets
data.dir = switch(.type,
                  biolog = file.path("biologicalData", input.data.dir),
                  synthet = file.path("syntheticData", input.data.dir))

if(.type == "synthet"){
  
  pattern = "bestWHlist"

  runs = .stats.runs
  sub.dir.files.tmp = sub.dir.files
  sub.dir.figures.tmp = sub.dir.figures
  sub.dir.RData.tmp = sub.dir.RData
  
  # run <- 1
  # SWITCH: over the NYC trip, we will run 1:43 for sCCA and 44:86 for NNMF.
  #   Then, on the morning of Tuesday, 25 June, we swap these indices. This
  #   ensures that both methods had read ability on the main data files.
  # for(run in 94:100){
  # for(run in 87:93){
  # for(run in 1:43){
  # for(run in 44:86){
  for(run in 1:runs){
    
    cat("Current stats run:", run, "\n")
    
    sub.dir.files = file.path(sub.dir.files.tmp, paste0("run", run))
    sub.dir.figures = file.path(sub.dir.figures.tmp, paste0("run", run))
    sub.dir.RData = file.path(sub.dir.RData.tmp, paste0("run", run))
    
    results = list.files(sub.dir.RData, pattern)
    NMF_result = list()
    
    # i <- 1
    for(i in 1:length(results)){
      # load bestWHlist 
      load(file.path(sub.dir.RData, results[i]))
      cat("Current data set: ", paste(bestWHlist$datasets, collapse = "\t"), "\n")
      
      # get desired number of features to be identified
      desired <- getDesiredSynthet(
        bestWHlist, input.data.dir, stats.run = paste0("run", run)
      )
      cat("Number of desired features: ", as.numeric(desired), "\n")

      fname.merge = Reduce(intersect, strsplit(bestWHlist$datasets, split = "synth"))
      current_res <- extractAndSave(
        bestWHlist, desired, fnames = bestWHlist$datasets, fname.add = fname.merge
      )

      NMF_result = c(NMF_result, current_res)

    }
    
    save(NMF_result, file = file.path(sub.dir.RData, paste0("NMF_result.RData")))
    
    drawROC(NMF_result, method = "NMF", input.data.dir, stats.run = paste0("run", run))
    cat("Finished drawing ROC curves.\n")
    

    cat("Calculating AUCs, ACC, F1-score, Cohen's Kappa and MCC ...\n")
    auc_arr = calculateAUC(NMF_result, method = "NMF", input.data.dir, 
                           stats.run = paste0("run", run))
    
    acc_arr = calculateACC(NMF_result, method = "NMF", input.data.dir,
                           stats.run = paste0("run", run))
    
    f1_arr = calculateACC(NMF_result, method = "NMF", input.data.dir,
                         stats.run = paste0("run", run), measure = "f1")
    
    kappa_arr = calculateACC(NMF_result, method = "NMF", input.data.dir,
                          stats.run = paste0("run", run), measure = "kappa")
    
    mcc_arr = calculateACC(NMF_result, method = "NMF", input.data.dir,
                          stats.run = paste0("run", run), measure = "mcc")

  }
  
  # store results of stats runs in one list
  NMF_result_list = vector("list", runs) 
  
  # for(run in 94:100){
  # for(run in 87:93){
  # for(run in 1:43){
  # for(run in 44:86){
  for(run in 1:runs){
    
    sub.dir.RData = file.path(sub.dir.RData.tmp, paste0("run", run))
    load(file.path(sub.dir.RData, "NMF_result.RData"))
    NMF_result_list[[run]] = NMF_result
    
  }
  
  saveRDS(NMF_result_list, file = file.path(sub.dir.RData.tmp, "NMF_result_list.RData"))
  
} else {
  
  pattern = "bestWHlist.RData"
  results = list.files(sub.dir.RData, pattern)
  NMF_result = list()
  
  for(i in 1:length(results)){
    # load bestWHlist 
    load(file.path(sub.dir.RData, results[i]))
    cat("Current data set: ", paste(bestWHlist$datasets, collapse = "\t"), "\n")
    
    # get desired number of features to be identified
    desired = getDesiredBiolog(bestWHlist) 
    cat("Number of desired features: ", as.numeric(desired), "\n")
    
    fname.merge = ""
    current_res = extractAndSave(bestWHlist, desired, fnames = bestWHlist$datasets,
                                 fname.add = fname.merge)
    
    NMF_result = c(NMF_result, current_res)
    
  }
  save(NMF_result, file = file.path(sub.dir.RData, paste0("NMF_result.RData")))
  
  # do cross validation 
  if(.do.CV == TRUE) {
    
    source(file.path(.src.dir, "crossValidation.R"))
    results_train = list.files(sub.dir.RData, "bestWHlist_train")
    
    for(s in 1:.subsets.val){
      
      # bestWHlist
      load(file.path(sub.dir.RData, results_train[s]))
      training = sub("bestWHlist_train", "", sub(".RData", "", results_train[s]))
      
      # get desired number of features to be identified
      desired = getDesiredBiolog(bestWHlist)
      cat("Training set", training, "\n",
          "Get ", unlist(desired), " top ranked weights in", bestWHlist$datasets, "\n")
      
      current_res = extractAndSave(bestWHlist, desired, fnames = bestWHlist$datasets,
                                   fname.add2 = paste0("_train", training))

      }
    
    evaluatePerformance("NMF", subsets = .subsets.val)
    
  }
}

cat("Finished calculation of post-processing of NMF.\n")

sink()
close(log.con)

  
  