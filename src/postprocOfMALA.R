#######################################################################
##### file: postprocOfMALA.R                                      #####
##### input: logic formulas for sample classification .csv        #####
##### output: plane list of cancidate genes                       #####
##### packages:                                                   #####
##### author: B. Pucher                                           #####
##### date created: 10/06/2016                                    #####
##### last change: 29/11/2017                                     #####
#######################################################################

rm(list = ls())

#install.packages("abind")
library(abind)

#######################################################################
#####                                                             #####
##### FUNCTIONS                                                   #####
#####                                                             #####
#######################################################################

#####------------------------------------------------------------------
# Extract MALA features from formulas and save results to .txt and .RData
#####------------------------------------------------------------------
extractAndSave = function(mala.out.dir, desired, fnames, fname.add = "",
                          fname.add2 = ""){
  
  MALA_res = extractMALAFeatures(mala.out.dir, desired)
  
  invisible(mapply(function(x, n) {
    write.table(x, file.path(sub.dir.files, paste0(n, ".txt")),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
  }, x = MALA_res, n = paste0("MALA_", fnames, fname.add2)))
  cat("Saved features selected by MALA to .txt-files. \n")
  
  saveRDS(MALA_res, file = file.path(sub.dir.RData, paste0("MALA_result", fname.add,
                                                           fname.add2, ".RData")))
  cat("Saved features selected by MALA to .RData-files. \n\n")
  
  return(MALA_res)
  
}

#####------------------------------------------------------------------
# Calculate and save global performance statistics manually. 
#####------------------------------------------------------------------
calculateGlobalStats = function(mala.out.path, set){
  
  pop.var = function(x) var(x) * (length(x)-1) / length(x)
  pop.sd = function(x) sqrt(pop.var(x))
  
  # read in stats file for all splits
  stats.files = list.files(mala.out.path, pattern = paste0(set, ".stats"), full.names = TRUE)
  
  stats.df.na = lapply(stats.files, read.csv, skip = 1, stringsAsFactors = FALSE, row.names = 1)
  stats.df = lapply(stats.df.na, function(x) Filter(function(df) !all(is.na(df)) , x))
  stats.arr = abind(stats.df, along = 3)
  
  # calculate global statistics
  stats.mean = apply(stats.arr, MARGIN = c(1,2), function(x) round(mean(x), 2))
  stats.var = apply(stats.arr, MARGIN = c(1,2), function(x) round(pop.var(x), 2))
  stats.sd = apply(stats.arr, MARGIN = c(1,2), function(x) round(pop.sd(x), 2))
  
  # write global statistics to file
  global.stats.file = file(file.path(mala.out.path, paste0(set, "globalstats_man.csv")), 
                           open = "w")
  writeLines(con = global.stats.file, text = "MEANS:")
  write.csv(stats.mean, file = global.stats.file, quote = FALSE)
  writeLines(con = global.stats.file, text = "VARIANCES:")
  write.csv(stats.var, file = global.stats.file, quote = FALSE)
  writeLines(con = global.stats.file, text = "STANDARD DEVIATIONS:")
  write.csv(stats.sd, file = global.stats.file, quote = FALSE)
  
  close(global.stats.file)
  
}

#####------------------------------------------------------------------
# Create overview of global stats for parameter sets.
#####------------------------------------------------------------------
globalStatsOverview = function(mala.out.dir, set, use.man = TRUE, 
                               out.path = sub.dir.files.tmp,
                               measure = "Percentage_of_correct_classified_elements:"){

  # read in stats file for current parameter set
  file.pattern = ifelse(use.man, yes = paste0(set, "globalstats_man"),
                        no = paste0(set, "globalstats.csv"))
  
  stats.files = sapply(file.path(out.path, paste0("run", 1:.stats.runs), 
                                 mala.out.dir), list.files, 
                       pattern = file.pattern, full.names = TRUE)
  
  stats.df = lapply(stats.files, read.csv, skip = 1, stringsAsFactors = FALSE,
                    row.names = 1, nrows = 7)
  names(stats.df) = paste0("run", 1:.stats.runs)
    
  stats.line = do.call(rbind, lapply(stats.df, function(df) {
    df[measure, ] }))
  
  return(list(stats.line))

}

#######################################################################
#####                                                             #####
##### MAIN SECTION                                                #####
#####                                                             #####
#######################################################################
current.run = .current.MALA
input.data.dir = switch(.type,
                        biolog = .current.biolog,
                        synthet = .current.synthet)
sub.dir.name = file.path(paste0(.type, "MALA"), current.run)
source(file.path(.src.dir, "setSubDirPath.R"))
source(file.path(.src.dir, "extractMALAFeatures.R"))
source(file.path(.src.dir, "drawROC.R"))
source(file.path(.src.dir, "assessPerformance.R"))

log.con = file(file.path(sub.dir.files, "postprocOfMALALog.txt"))
sink(file = log.con, type = "output")
flush(log.con)

cat("Start postprocessing of MALA ... \n")
cat("Assuming datasets from run:", input.data.dir, "\n")

# load input datasets for cross validation
data.dir = switch(.type,
                  biolog = file.path("biologicalData", input.data.dir),
                  synthet = file.path("syntheticData", input.data.dir))

if(.type == "synthet"){

  runs = .stats.runs
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
    
    datasets.all = list.files(file.path(.run.dir, "syntheticData", input.data.dir, 
                                        "RData", paste0("run", run)))
    datasets = datasets.all[grep("_Diff|train|test|PWs", datasets.all, invert = TRUE)] 
    
    # get n_groups, p_de, p_dm, c_de, c_dm from datasets names
    source(file.path(.src.dir, "helpParameter.R")) 
    
    MALA_result = vector("list", length(n_groups)*length(p_de)*length(c_de)*2)
    r = 1
    
    for(k in 1:length(n_groups)){
      for(i in 1:length(p_de)){
        for(j in 1:length(c_de)){
          
          fname_ge = paste("GeneExp_synth_", n_groups[k], "gr_", p_de[i], "p_", c_de[j], sep = "")
          fname_met = paste("Methyl_synth_",  n_groups[k], "gr_", p_dm[i], "p_", c_dm[j], sep = "")
          data.names = c(fname_ge, fname_met)
          
          params = paste0(letters[k], letters[i], letters[j])
          
          mala.out.dir = paste(.type, params, "out", sep = "_")
          
          desired = sum(unlist(getDesiredSynthet(input.data.dir, run = run)))
          
          # get names of features used in logic formulas
          MALA_features = extractMALAFeatures(mala.out.dir, desired)
          MALA_genes = MALA_features$GeneExp
          MALA_met = MALA_features$Methyl
          MALA_result[[r]] <- MALA_genes; names(MALA_result)[r] <- fname_ge
          MALA_result[[r+1]] <- MALA_met; names(MALA_result)[r+1] <- fname_met
          r = r+2
          
          write.table(MALA_genes, file.path(sub.dir.files, paste0("MALA_", fname_ge, ".txt")),
                      quote = FALSE, row.names = FALSE, col.names = FALSE)
          write.table(MALA_met, file.path(sub.dir.files, paste0("MALA_", fname_met, ".txt")),
                      quote = FALSE, row.names = FALSE, col.names = FALSE)
          
          cat("Saved features detected by sCCA to files.\n\n")
 
          # calculate global statistics for current mala.out.dir
          mala.out.path = file.path(sub.dir.files, mala.out.dir)
          
          calculateGlobalStats(mala.out.path, set = "TRAIN")
          calculateGlobalStats(mala.out.path, set = "TEST")
          
        }
      }
    }
    
    cat("Finished calculation of MALA.\n")
    
    save(MALA_result, file = file.path(sub.dir.RData, "MALA_result.RData"))
    
    drawROC(MALA_result, method = "MALA", input.data.dir, stats.run = paste0("run", run))
    cat("Finished drawing ROC curves.\n")
    
    cat("Calculating AUCs, ACC, F1-score, Cohen's Kappa and MCC ...\n")
    auc_arr = calculateAUC(MALA_result, method = "MALA", input.data.dir, 
                           stats.run = paste0("run", run))

    acc_arr = calculateACC(MALA_result, method = "MALA", input.data.dir, 
                           stats.run = paste0("run", run))
    
    f1_arr = calculateACC(MALA_result, method = "MALA", input.data.dir,
                         stats.run = paste0("run", run), measure = "f1")
    
    kappa_arr = calculateACC(MALA_result, method = "MALA", input.data.dir,
                         stats.run = paste0("run", run), measure = "kappa")
    
    mcc_arr = calculateACC(MALA_result, method = "MALA", input.data.dir,
                         stats.run = paste0("run", run), measure = "mcc")
    
  }
  
  # store results of stats runs in one list
  MALA_result_list = vector("list", runs) 
  for(run in 1:runs){
    
    sub.dir.RData = file.path(sub.dir.RData.tmp, paste0("run", run))
    load(file.path(sub.dir.RData, "MALA_result.RData"))
    MALA_result_list[[run]] = MALA_result
  }
  
  saveRDS(MALA_result_list, file = file.path(sub.dir.RData.tmp, "MALA_result_list.RData"))
  
  sub.dir.files = sub.dir.files.tmp
  sub.dir.figures = sub.dir.figures.tmp
  sub.dir.RData = sub.dir.RData.tmp
  
  # Create overview of global stats for simulation parameter sets
  overview.TRAIN = list()
  overview.TEST = list()
    
  for(k in 1:length(n_groups)){
    for(i in 1:length(p_de)){
      for(j in 1:length(c_de)){
        
        params = paste0(letters[k], letters[i], letters[j])
        mala.out.dir = paste(.type, params, "out", sep = "_")
        
        # overview of stats on TRAIN sets
        current.ov.TRAIN = globalStatsOverview(mala.out.dir, set = "TRAIN", 
                                               out.path = sub.dir.files)
        names(current.ov.TRAIN) = mala.out.dir
        overview.TRAIN = c(overview.TRAIN, current.ov.TRAIN)

        
        # overview of stats on TEST sets
        current.ov.TEST = globalStatsOverview(mala.out.dir, set = "TEST", 
                                               out.path = sub.dir.files)
        names(current.ov.TEST) = mala.out.dir
        overview.TEST = c(overview.TEST, current.ov.TEST)
        
      }
    }
  }
  
  # merge overviews into data frames
  overview.TRAIN.df = do.call(cbind, overview.TRAIN)
  overview.TRAIN.df = Filter(function(x) !all(x == 0), overview.TRAIN.df)
  colnames(overview.TRAIN.df) = sub("(.*)\\.", "", colnames(overview.TRAIN.df))
    
  overview.TEST.df = do.call(cbind, overview.TEST)
  overview.TEST.df = Filter(function(x) !all(x == 0), overview.TEST.df)
  colnames(overview.TEST.df) = sub("(.*)\\.", "", colnames(overview.TEST.df))
  
  # write overview to files
  ov.TRAIN.file = file(file.path(sub.dir.files, "TRAINglobalstats_overview.csv"), open = "w")
  writeLines(con = ov.TRAIN.file, text = paste0(",", names(overview.TRAIN), ",", collapse = ","))
  write.csv(overview.TRAIN.df, file = ov.TRAIN.file, quote = FALSE)
  close(ov.TRAIN.file)
  
  ov.TEST.file = file(file.path(sub.dir.files, "TESTglobalstats_overview.csv"), open = "w")
  writeLines(con = ov.TEST.file, text = paste0(",", names(overview.TEST), ",", collapse = ","))
  write.csv(overview.TEST.df, file = ov.TEST.file, quote = FALSE)
  close(ov.TEST.file)
  
        
} else {
  
  mala.out.dir = "biolog_out"

  datasets.all = list.files(file.path(.run.dir, "biologicalData", input.data.dir, "RData"))
  datasets = datasets.all[grep("train|test", datasets.all, invert = TRUE)]
  data.names = gsub(".RData", "", datasets)
  data.paths = as.list(file.path(.run.dir, "biologicalData", input.data.dir, "RData", datasets))
  data.raw = lapply(data.paths, readRDS)
  
  desired = sum(unlist(getDesiredBiolog(data.raw, data.names)))
  cat("Number of desired features: ", as.numeric(desired), "\n")
  
  MALA_result = extractAndSave(mala.out.dir, desired, fnames = data.names)

  cat("Number of selected features in each data type:", 
      mapply(paste, names(MALA_result), lapply(MALA_result, length)), "\n")
  cat("Number of selected features in total:", length(unlist(MALA_result)), "\n")
  
  # post-processing of cross validation 
  if(.do.CV == TRUE){
    
    source(file.path(.src.dir, "crossValidation.R"))
    
    for(s in 1:.subsets.val){
    
      mala.out.dir = paste0("biolog_train", s, "_out")
      
      cat("Training set", s, "\n",
          "Get ", unlist(desired), " features in up to 100 MALA runs. \n")
      
      current_res = extractAndSave(mala.out.dir, desired, fnames = data.names,
                                   fname.add2 = paste0("_train", s))
      
      cat("Number of selected features in each data type:", 
          mapply(paste, names(current_res), lapply(current_res, length)), "\n")
      cat("Number of selected features in total:", length(unlist(current_res)), "\n")
    }
    
    evaluatePerformance("MALA", subsets = .subsets.val)
  }
}
    
cat("Finished calculation of post-processing of MALA.\n")    
cat("*** MALA DONE ***")


sink()
close(log.con)
