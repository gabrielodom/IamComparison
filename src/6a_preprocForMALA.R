#######################################################################
##### file: preprocForMALA.R                                      #####
##### input: /results/run/*biological or synthetic datasets*      #####
#####     in seperate txt files                                   #####
##### output: concatenated datasets in DMBCSV format expected     #####
#####     by MALA                                                 #####
##### packages: --                                                #####
##### author: B. Pucher                                           #####
##### date created: 16/06/2015                                    #####
##### last change:  09/11/2017                                    #####
#######################################################################

rm(list = ls())

#####------------------------------------------------------------------
# TCGA barcode structure labeling samples
#####------------------------------------------------------------------
barcode.parts = c("Project", "TSS", "Participant", "Sample", "Analyte",
                  "Plate", "Center")

#####------------------------------------------------------------------
# tissue type code from TCGA Wiki at National Cancer Institute
#####------------------------------------------------------------------
tumor = 1
normal = 11
metastatic = 6
control = 20

#######################################################################
#####                                                             #####
##### FUNCTIONS                                                   #####
#####                                                             #####
#######################################################################

#####------------------------------------------------------------------
# rename rows (Features) of biologic dataset
#####------------------------------------------------------------------
renameSaFt = function(dataset){
  sa.names = colnames(dataset)
  sa.type = sub("_.*", "", sa.names)
  analyte = unique(sub(".*_", "", sa.names))
  switch(analyte, 
         R = {
           add = "ge."
         },
         D = {
           add = "met."
         }
  )
  rownames(dataset) = paste(add, rownames(dataset), sep = "")
  return(dataset)
}


#####------------------------------------------------------------------
# create csv file in DMB-format
#####------------------------------------------------------------------
createDMBCSV = function(conc.data, fname){
  class.row = switch(.type,
                     biolog = c("class", "", sub("_.*", "", colnames(conc.data))),
                     synthet = c("class", "", sub(".*_", "", colnames(conc.data))))
  
  data.df = cbind(gene = rownames(conc.data), 
                  type = rep("NUM", nrow(conc.data)), conc.data)
  
  csv.con = file(file.path(sub.dir.files, fname), open = "w+")
  write.csv(rbind(class.row, data.df), csv.con, quote = F,
            row.names = F)
  close(csv.con)
  cat("Dataset standardized and saved in DMBCSV format to file", 
      fname, "\n")
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

log.con = file(file.path(sub.dir.files, "preprocForMALALog.txt"))
sink(file = log.con, type = "output")
flush(log.con)

cat("Start preprocessing", .type, "datasets for MALA ... \n")
cat("Datasets created in run:", input.data.dir, "\n")

# load datasets
data.dir = switch(.type,
                  biolog = file.path("biologicalData", input.data.dir),
                  synthet = file.path("syntheticData", input.data.dir))
datasets.all = list.files(file.path(.run.dir, data.dir, "RData"))

if(.type == "synthet") { 
  
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
    
    datasets.all = list.files(file.path(.run.dir, data.dir, "RData", paste0("run", run)))
    datasets = datasets.all[grep("_Diff|train|test|PWs", datasets.all, invert = TRUE)] 
    
    # get n_groups, p_de, p_dm, c_de, c_dm from datasets names
    source(file.path(.src.dir, "helpParameter.R"))   
    
    for(k in 1:length(n_groups)){
      for(i in 1:length(p_de)){
        for(j in 1:length(c_de)){
          
          fname_ge = paste("GeneExp_synth_", n_groups[k], "gr_", p_de[i], "p_", c_de[j], sep = "")
          fname_met = paste("Methyl_synth_",  n_groups[k], "gr_", p_dm[i], "p_", c_dm[j], sep = "")
          data.names = c(fname_ge, fname_met)
          data.paths = as.list(file.path(.run.dir, data.dir, "RData", paste0("run", run), 
                                         paste0(data.names, ".RData")))
          data.raw = lapply(data.paths, readRDS)
          names(data.raw) = data.names
          
          # standardize the rows (features) of the data sets
          data.st = lapply(lapply(lapply(data.raw, t), scale), t)
          
          data = data.st
          
          # concatenate datasets in such a way that the values of all features 
          # measured on one sample show up in one column (one col per sample)
          data.whole = do.call(rbind, data)
          
          # create data frame with adequat (DMB) header and save to csv-file 
          #params = paste(n_groups[k], unique(c(p_de[i], p_dm[i])), unique(c(c_de[j], c_dm[j])), sep = "_", collapse = "-")
          params = paste0(letters[k], letters[i], letters[j])
          file.name = paste(.type, params, "DMB.csv", sep = "_")
          createDMBCSV(data.whole, file.name)
        }
      }
    }
  }
  
} else {
  
  datasets = datasets.all[grep("train|test", datasets.all, invert = TRUE)]
  data.names = gsub(".RData", "", datasets)
  data.paths = as.list(file.path(.run.dir, data.dir, "RData", datasets))
  data.raw = lapply(data.paths, readRDS)
  names(data.raw) = data.names
  
  # standardize the rows (features) of the data sets
  data.st = lapply(lapply(lapply(data.raw, t), scale), t)
  
  # rename rows (features) to contain the data type also 
  data = lapply(data.st, renameSaFt)

  # concatenate datasets in such a way that the values of all features 
  # measured on one sample show up in one column (one col per sample)
  data.whole = do.call(rbind, data)
  
  # create data frame with adequat (DMB) header and save to csv-file 
  file.name = paste0(.type, "_DMB.csv")
  createDMBCSV(data.whole, file.name)
  
  # for cross validation split dataset into training- and tets set
  if(.do.CV == TRUE){
    
    source(file.path(.src.dir, "crossValidation.R"))
    
    # preprocess training data
    for(s in 1:.subsets.val){
      
      data.path.train = file.path(.run.dir, data.dir, "RData", paste0("data_train", s, ".RData"))
      data.train.raw = readRDS(data.path.train)
      
      # standardize the rows (features) of the data sets
      data.train.st = lapply(lapply(lapply(data.train.raw, t), scale), t)
      
      # rename rows (features) to contain the data type also 
      data.train = lapply(data.train.st, renameSaFt)

      # concatenate datasets in such a way that the values of all features 
      # measured on one sample show up in one column (one col per sample)
      data.train.whole = do.call(rbind, data.train)
      
      # create data frame with adequat (DMB) header and save to csv-file 
      file.name = paste0(.type, "_train", s, "_DMB.csv")
      createDMBCSV(data.train.whole, file.name)
    }
  }
}

cat("Preprocessing of", .type, "datasets finished. \n")

sink()
close(log.con)


