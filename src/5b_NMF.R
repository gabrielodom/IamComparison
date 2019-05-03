#######################################################################
##### file: 5b_NMF.R                                              #####
##### input: *.csv files adequat for input to NMF from preproc.   #####
##### output: best factor matrices W and H to rectonstuct data    #####
##### packages: --                                                #####
##### author: B. Pucher                                           #####  
##### date created: 13/04/2015                                    #####
##### last change:  09/11/2017                                    #####
#######################################################################
rm(list = ls())

#######################################################################
#####                                                             #####
##### FUNCTIONS                                                   #####   
#####                                                             #####
#######################################################################

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
source(file.path(.src.dir, "multiNMF_mm.R"))
source(file.path(.src.dir, "multiNMF_residue.R"))
source(file.path(.src.dir, "multiNMF_residue_comodule.R"))

log.con = file(file.path(sub.dir.files, "NMFLog.txt"))
sink(file = log.con, type = "output")
flush(log.con)

cat("Assuming datasets created in run:", input.data.dir, "\n")

if(.type == "synthet"){
  
  nloop = .nloop
  maxiter = .maxiter
  
  # run NMF with number of basis vectors K 
  K = .good.k
  
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
    
    dir.create(sub.dir.files, showWarnings = FALSE)
    dir.create(sub.dir.figures, showWarnings = FALSE)
    dir.create(sub.dir.RData, showWarnings = FALSE)
    
    # get number of nn datasets in synthetic data dir)
    datasets_nn = sort(list.files(sub.dir.RData, "data_nn"))
    
    run.seed = stats.seed[run]
    cat("Seed used:", run.seed, "\n")
    set.seed(run.seed)
    
    for(i in 1:length(datasets_nn)){
      load(file.path(sub.dir.RData, datasets_nn[i]))
      
      cat("Datasets: ", paste(names(data), collapse = "\t"), "\n")
      cat("calcualte NMF with K = ", K, " nloop = ", nloop, " and maxiter = ", maxiter, "\n")
      cat(paste("---------------", Sys.time(), "---------------"), "\n")
      bestWHlist_temp = multiNMF_residue_comodule(data, K, nloop, maxiter)
      bestWHlist = c(bestWHlist_temp, datasets = list(names(data)))
      save(bestWHlist,
           file = file.path(sub.dir.RData, paste0("bestWHlist", unlist(strsplit(datasets_nn[i], split = "data_nn"))[2])))
      cat(paste("---------------", Sys.time(), "---------------"), "\n")
    
    }
  }
  
} else {
  
  # get number of nn datasets in biologic data dir)
  datasets_nn_all = list.files(sub.dir.RData, "data_nn")
  datasets_nn_train = datasets_nn_all[grep("train", datasets_nn_all)]
  datasets_nn = datasets_nn_all[grep("train", datasets_nn_all, invert = TRUE)]
  
  for(i in 1:length(datasets_nn)){
    load(file.path(sub.dir.RData, datasets_nn[i]))
    
    # run NMF with number of basis vectors K 
    nloop = .nloop
    maxiter = .maxiter
    K = .good.k
    
    set.seed(.biolog.NMF.seed)
    cat("Datasets: ", paste(names(data), collapse = "\t"), "\n")
    cat("Number of samples in data set: ", unique(unlist(lapply(data, nrow))), "\n")
    cat("calcualte NMF with K = ", K, " nloop = ", nloop, " and maxiter = ", maxiter, "\n")
    cat(paste("---------------", Sys.time(), "---------------"), "\n")
    bestWHlist_temp = multiNMF_residue_comodule(data, K, nloop, maxiter)
    bestWHlist = c(bestWHlist_temp, datasets = list(names(data)))
    save(bestWHlist,
         file = file.path(sub.dir.RData, paste0("bestWHlist", unlist(strsplit(datasets_nn[i], split = "data_nn"))[2])))
    cat(paste("---------------", Sys.time(), "---------------"), "\n")
    
  }; rm(data)
  
  if(.do.CV == TRUE){
    
    # run NMF with number of basis vectors K
    nloop = .nloop
    maxiter = .maxiter
    K = .good.k
    
    for(i in 1:length(datasets_nn_train)){
      load(file.path(sub.dir.RData, datasets_nn_train[i])) #data
      # get current training set number from file name
      training = sub("data_nn_train", "", sub(".RData", "", datasets_nn_train[i]))
      cat("---Split number---", training, "\n")
      cat("Datasets: ", paste(names(data), collapse = "\t"), "\n")
      cat("Number of samples in training set: ", unique(unlist(lapply(data, nrow))), "\n")
      cat("calcualte NMF with K = ", K, " nloop = ", nloop, " and maxiter = ", maxiter, "\n")
      cat(paste("-----", Sys.time(), "-----"), "\n")
      bestWHlist_temp = multiNMF_residue_comodule(data, K, nloop, maxiter)
      bestWHlist = c(bestWHlist_temp, datasets = list(names(data)))
      
      save(bestWHlist, file = file.path(sub.dir.RData, paste0("bestWHlist",
                                                              sub("data", "", sub("_nn", "", datasets_nn_train[i])))))   
      cat(paste("-----", Sys.time(), "-----"), "\n")
      
    }
  } 
  
}


cat("Finished calculation of NMF.\n")  
  
sink()
close(log.con)


