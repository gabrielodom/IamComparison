#######################################################################
##### Download TCGA data with TCGA-Assembler Version 2            #####
#######################################################################

#install.packages("httr")
#install.packages("HGNChelper")
#install.packages("RCurl")
#install.packages("rjson")
#install.packages("stringr")

dest = .data.raw
  
# Load functions
source(file.path(.TCGA.assembler.dir, "Module_A.R"))
source(file.path(.TCGA.assembler.dir, "Module_B.R"))

#' set data saving path
sPath0 = dest
sPath1 = file.path(dest, "rawData")
sPath2 = file.path(dest, "prbasicData")
sPath3 = file.path(dest, "pradvData")

#' choose a cancer type
sCancer = .sCancer

#######################################################################
##### Part 1: Downloading Data of different platforms using Module A 
#######################################################################

# Download DNA methylation 450 data
path_methylation_450 <-
  DownloadMethylationData(cancerType = sCancer,
                          assayPlatform = "methylation_450",
                          saveFolderName = sPath1)

# Download DNA methylation 27 data
path_methylation_27 <-
  DownloadMethylationData(cancerType = sCancer,
                          assayPlatform = "methylation_27",
                          saveFolderName = sPath1)


# Download gene expression data
path_geneExp <-
  DownloadRNASeqData(cancerType = sCancer,
                     assayPlatform = "gene.normalized_RNAseq",
                     saveFolderName = sPath1)

#######################################################################
#### Part 2: Perform basic processing of downloaded data using Module B 
#######################################################################

# Process DNA methylation 450 data
list_methylation_450 <-
  ProcessMethylation450Data(inputFilePath = path_methylation_450[1],
                            outputFileName = paste(sCancer, "methylation_450",
                                                   sep = "__"),
                            outputFileFolder = sPath2)

# Process DNA methylation 27 data
list_methylation_27 <-
  ProcessMethylation27Data(inputFilePath = path_methylation_27[1],
                            outputFileName = paste(sCancer, "methylation_27",
                                                   sep = "__"),
                            outputFileFolder = sPath2)

# Merge DNA methylation 450 and 27 data
list_methylation_merged =
  MergeMethylationData(input1 = list_methylation_27, 
                       input2 = list_methylation_450,
                       outputFileName = paste(sCancer, "methylatioin27_450_merged",
                                              sep = "__"),
                       outputFileFolder = sPath2)

# Process gene expression data
list_geneExp <-
  ProcessRNASeqData(inputFilePath = path_geneExp[1],
                    outputFileName = paste(sCancer, "geneExp", sep = "__"),
                    dataType = "geneExp",
                    outputFileFolder = sPath2)

#######################################################################
#### Part 3: Create inputDataList with downloaded data.  
#######################################################################

# 1st step: Load multiplatform data
l_methylation     <- list(Des  = list_methylation_merged$Des,
                          Data = list_methylation_merged$Data,
                          dataType = "Methylation")
l_geneExp         <- list(Des  = list_geneExp$Des,
                          Data = list_geneExp$Data,
                          dataType = "GeneExp")

# 2nd step: Put multiplatform data in a vector of list objects
inputDataList      <- vector("list", 2)
inputDataList[[1]] <- l_geneExp
inputDataList[[2]] <- l_methylation

save(inputDataList, file = file.path(sPath0, paste0("tcga_", tolower(sCancer), 
                                                    ".RData")))




