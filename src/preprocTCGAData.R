#######################################################################
##### file: preprocTCGAData.R                                     #####
##### input: TCGA data .dataset.list stored at data/raw/          #####
##### output: results/*CurrentRun*/*SubDir*/RData/subsets.RData   #####
##### packages: --                                                #####
##### author: B. Pucher                                           #####
##### date created: 13/07/2015                                    #####
##### last change:  09/11/2017                                    #####
#######################################################################
rm(list=ls())

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
# remove all data types not specified from inputDataList
#####------------------------------------------------------------------
keepDataTypes = function(dataset.list, types.to.keep){
  keep = logical(length(dataset.list)) #initially all FALSE
  for(i in 1:length(dataset.list)){
    if(dataset.list[[i]]$dataType %in% types.to.keep){
      keep[i] = TRUE
    } 
  }
  
  data.list = dataset.list
  data.list[!keep] <- NULL
  return(data.list)
}

#####------------------------------------------------------------------
# extract data type of each data set
#####------------------------------------------------------------------
getDataType = function(dataset){
  dataset$dataType
}

#####------------------------------------------------------------------
# print out number of features and number of samples in a dataset
#####------------------------------------------------------------------
printDimension = function(dataset){
  cat(dataset$dataType, "Features:", dim(dataset$Data)[1],
        "Samples:", dim(dataset$Data)[2], "\n")
}

#######################################################################
##### Feature Manipulations                                       #####   
#######################################################################

#####------------------------------------------------------------------
# extract gene names (gene symbols) of a data set
#####------------------------------------------------------------------
getGeneSymbols = function(dataset){
  dataset$Des[,"GeneSymbol"]
}

#####------------------------------------------------------------------
# extract sample names of a data set
#####------------------------------------------------------------------
getSampleNames = function(dataset){
  dimnames(dataset$Data)[[2]]
}

#####------------------------------------------------------------------
# split sample names of TCGA datasets into 7 blocks
# example sample name: "TCGA-A8-A07B-01A-11R-A007Z-07"
#####------------------------------------------------------------------
splitNames = function(names, sep = "-"){
  blocks = matrix(unlist(strsplit(names, split = sep)), 
                  nrow = length(names), byrow = TRUE,
                  dimnames = list(NULL, barcode.parts)) 
  
  return(blocks)
}

#####------------------------------------------------------------------
# split the 4th block of sample names to access tissue type
#####------------------------------------------------------------------
sampleType = function(name.blocks){
  type = matrix(unlist(strsplit(name.blocks[,"Sample"], split = "")),
                nrow = nrow(name.blocks), byrow = TRUE)
  type.fc = factor(as.numeric(paste(type[,1], type[,2], sep = "")),
                   levels = c(tumor, metastatic, normal, control),
                   labels = c("tumor", "metastatic", "normal", "control"))
  return(type.fc)
}

#####------------------------------------------------------------------
# table with tumor and normal samples per data type
#####------------------------------------------------------------------
getSampleTypes = function(dataset){
  # names of genes and samples in the data set
  sample.names = getSampleNames(dataset)
  gene.symbols = getGeneSymbols(dataset)
  
  # split sample names and get the tissue type of each sample 
  sample.name.blocks = splitNames(sample.names)
  sample.types = sampleType(sample.name.blocks)
  return(sample.types) 
}

typeTable = function(sample.types){
  return(table(sample.types))
}

#####------------------------------------------------------------------
# remove all elements of $Data that do not belong to the specified
# tissue type indicated by the index-subset
#####------------------------------------------------------------------
dataSubset = function(dataset){
  sample.types = getSampleTypes(dataset)
  tu.idx = which(as.character(sample.types) == "tumor") 
  no.idx = which(as.character(sample.types) == "normal")
  tu.subset = dataset
  no.subset = dataset
  tu.subset$Data = tu.subset$Data[,tu.idx]
  tu.subset$dataType = paste(dataset$dataType, "_tu", sep = "")
  no.subset$Data = no.subset$Data[,no.idx]
  no.subset$dataType = paste(dataset$dataType, "_no", sep = "")
  return(list(tu.subset, no.subset))
}

#######################################################################
#####                                                             #####
##### MAIN SECTION                                                #####
#####                                                             #####
#######################################################################
sub.dir.name = "preprocTCGAData-2"
source(file.path(.src.dir, "setSubDirPath.R"))
       
log.con = file(file.path(sub.dir.files, "preprocTCGALog.txt"))
sink(file = log.con, type = "output")
flush(log.con)

load(file.path(.data.raw, .dataset.list)) # inputDataList

# get rid of data sets not specified in parameters file 
data.list.raw = keepDataTypes(inputDataList, .data.to.integrate)
                            
cat("datasets raw: \n")
invisible(lapply(data.list.raw, printDimension))
rm(inputDataList)

data.symbol.complete = data.list.raw

# store lookup table of feature identifiers and gene symbols
invisible(lapply(data.symbol.complete, function(x) {
  saveRDS(x$Des, file = file.path(sub.dir.RData, paste0(x$dataType, "_Des.RData")))
}))

invisible(lapply(data.symbol.complete, function(x) {
  write.table(x$Des, file = file.path(sub.dir.RData, paste0(x$dataType, "_Des.tab")),
              quote = FALSE, row.names = FALSE, sep = "\t")
}))

# get type of each sample per data set
sample.types = lapply(data.symbol.complete, getSampleTypes)
names(sample.types) = unlist(lapply(data.symbol.complete, getDataType))

# list number of samples per tissue type in each data set
print(lapply(sample.types, typeTable))

# get one subset comprising tumor and one subset comprising normal 
# samples (if available) per data set
tu.no.subsets = lapply(data.symbol.complete, dataSubset)
cat("subsets of tumor and normal tissue samples: \n")
invisible(lapply(unlist(tu.no.subsets, recursive = F), printDimension))

subsets = unlist(tu.no.subsets, recursive = FALSE)
save(subsets, file = file.path(sub.dir.RData, "subsets.RData"))

# close the file connection and stop sinking
sink()
close(log.con)



