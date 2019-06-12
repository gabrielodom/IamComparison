#######################################################################
##### file: 5a_preprocForNMF.R                                    #####
##### input: biological or synthetic datasets located in          #####
#####        ./results/.current.run/biologicalData|syntheticData  #####
##### output: *.csv files adequat for input to NMF                #####
##### packages: --                                                #####
##### author: B. Pucher                                           #####
##### date created: 21/07/2015                                    #####
##### last change:  09/11/2017                                    #####
#######################################################################
rm(list = ls())

#######################################################################
#####                                                             #####
##### FUNCTIONS                                                   #####   
#####                                                             #####
#######################################################################

#####------------------------------------------------------------------
# calculate sum of squares of a matrix
#####------------------------------------------------------------------
sumOfSquares = function(matrix){
  return(sum(abs(matrix)^2))
}

#####------------------------------------------------------------------
# calculate the Frobenius norm of a matrix
#####------------------------------------------------------------------
frobeniusNorm = function(matrix){
  return(sqrt(sum(abs(matrix)^2)))
}

#####------------------------------------------------------------------
# scale data matrix to have Frobenius norm equal to a given number
#####------------------------------------------------------------------
adjustNorm = function(matrix, large.number = 5000){
  sosq = sumOfSquares(matrix)
  equal.sosq = large.number^2
  scale.factor = large.number^2/sosq
  matrix.sq.adj = matrix^2 * scale.factor
  matrix.adj = sqrt(matrix.sq.adj) * sign(matrix) 
  return(matrix.adj)
}

#####------------------------------------------------------------------
# fit data matrices to non-negativity constraints using method 
# suggested by Kim and Tidor (2003)
#####------------------------------------------------------------------
splitColumns = function(mtx){
  doubled.mtx = matrix(0, nrow = nrow(mtx), ncol = ncol(mtx)*2) 
  doubled.nms = character(length(dimnames(mtx)[[2]])*2)
  
  for(i in 1:ncol(mtx)){
    for(j in 1:nrow(mtx)){
      if(mtx[j,i] < 0){
        doubled.mtx[j, i*2] = abs(mtx[j,i])
      }
      else doubled.mtx[j,2*i-1] = mtx[j,i] 
    }
    
    doubled.nms[2*i-1] = colnames(mtx)[i]
    doubled.nms[2*i] = paste(doubled.nms[2*i-1], ".", sep = "")
  }
  
  dimnames(doubled.mtx) = list(rownames(mtx), doubled.nms)
  return(doubled.mtx)
}

#####------------------------------------------------------------------
# run preprocessing steps on dataset
#####------------------------------------------------------------------
doPreprocessing = function(data.raw, data.names){
  
  # standardize the columns (features) of the data matrices 
  # (= normalization according to Zhang (2012))
  data.st = lapply(data.raw, scale)
  cat("Features standardized. \n")
  
  # # scale matrices to have Frobenius Norm all equal to the same number
  data.fn = lapply(data.st, frobeniusNorm)
  cat("Frobenius norm of datasets .. \n", unlist(data.fn), "\n")
  
  data = lapply(data.st, adjustNorm, max(unlist(data.fn)))
  data.fn.eq = lapply(data, frobeniusNorm)
  cat("Frobenius norm after equalization .. \n", unlist(data.fn.eq), "\n")
  
  # fit data to non-negativity constraints 
  data.nn = lapply(data, splitColumns)
  names(data.nn) = data.names
  
  return(data.nn)
}

#####------------------------------------------------------------------
# save preprocessed dataset to .RData- and .csv-files 
#####------------------------------------------------------------------
savePreprocessed = function(data, fnames, save.csv = TRUE,
                            fname.add = "", fname.add2 = ""){
  if(save.csv == TRUE){
    mapply(function(x,n) {
      write.table(x, file.path(sub.dir.files, paste0(n, ".csv")),
                  quote = FALSE, row.names = TRUE, col.names = TRUE)
    }, x = data, n = paste0("nn_", fnames, fname.add2))
    cat("Non-negative dataset saved to csv files.\n")
  }
  
  save(data, file = file.path(sub.dir.RData, paste0("data_nn", fname.add, 
                                                    fname.add2, ".RData")))
  
  cat("Non-negative dataset saved to RData file. \n")
  
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

log.con = file(
  file.path(sub.dir.files, paste0("preprocForNMFLog_", Sys.Date(), ".txt"))
)
sink(file = log.con, type = "output")
flush(log.con)

cat("Start preprocessing datasets for analysis with NMF ..\n")
cat("Datasets created in run:", input.data.dir, "\n")

# load datasets
data.dir = switch(.type,
                  biolog = file.path("biologicalData", input.data.dir),
                  synthet = file.path("syntheticData", input.data.dir))

datasets.all = list.files(file.path(.run.dir, data.dir, "RData"))

if(.type == "synthet"){
  
  runs = .stats.runs
  sub.dir.files.tmp = sub.dir.files
  sub.dir.figures.tmp = sub.dir.figures
  sub.dir.RData.tmp = sub.dir.RData
  
  # run <- 1
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
    
    # k <- i <- j <- 1
    for(k in 1:length(n_groups)){
      for(i in 1:length(p_de)){
        for(j in 1:length(c_de)){
          
          fname_ge = paste("GeneExp_synth_", n_groups[k], "gr_", p_de[i], "p_", c_de[j], sep = "")
          fname_met = paste("Methyl_synth_",  n_groups[k], "gr_", p_dm[i], "p_", c_dm[j], sep = "")
          data.names = c(fname_ge, fname_met)
          data.paths = as.list(file.path(.run.dir, data.dir, "RData", paste0("run", run),
                                         paste0(data.names, ".RData")))
          data.raw.t = lapply(data.paths, readRDS)
          
          cat("Current datasets:", paste(data.names, collapse = "\t"), "\n")
          
          # transpose datasets 
          # NMF expects the features to be on the columns
          data.raw = lapply(data.raw.t, t)
          
          # standardize, euqalize Frobenius norm and fit to non-negativity
          data.nn = doPreprocessing(data.raw, data.names)
          
          # save datasets to *.RData and *.csv files
          fname.merge = Reduce(intersect, strsplit(data.names, split = "synth"))
          savePreprocessed(data.nn, save.csv = FALSE, fnames = data.names, fname.add = fname.merge)
          
        }
      }
    }
    
  }
  
} else {
  
  datasets = datasets.all[grep("train|test", datasets.all, invert = TRUE)] 
  data.names = gsub(".RData", "", datasets)
  data.paths = as.list(file.path(.run.dir, data.dir, "RData", datasets))
  data.raw.t = lapply(data.paths, readRDS)
  
  # transpose datasets 
  # NMF expects the features to be on the columns
  data.raw = lapply(data.raw.t, t)
  
  # standardize, euqalize Frobenius norm and fit to non-negativity
  data.nn = doPreprocessing(data.raw, data.names)
  
  # save datasets to *.RData and *.csv files
  savePreprocessed(data.nn, fnames = data.names)

  # for cross validation split dataset into training- and tets set
  if(.do.CV == TRUE) {
    
    source(file.path(.src.dir, "crossValidation.R"))
    
    # preprocess training data
    for(s in 1:.subsets.val){
      data.path.train = file.path(.run.dir, data.dir, "RData", paste0("data_train", s, ".RData"))
      data.train.t = readRDS(data.path.train)
      
      # transpose datasets
      data.train = lapply(data.train.t, t)
      
      # standardize, euqalize Frobenius norm and fit to non-negativity
      data.train.nn = doPreprocessing(data.train, names(data.train))
      
      # save datasets to *.RData and *.csv files
      savePreprocessed(data.train.nn, fnames = names(data.train), fname.add2 = paste0("_train", s))
      
    }
  }
}


cat("Finished preprocessing datasets for NMF.\n")

sink()
close(log.con)


