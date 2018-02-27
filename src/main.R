#######################################################################
##### Project: IntMethodCompPublication                           #####
##### main.R                                                      #####
#######################################################################
rm(list = ls(all.names = TRUE))

# Please set path to the Project folder 
.project.dir = "D:/Development/IntMethodCompPublication"

# Please enter the name of the list of datasets and 
# put the R-Object into .project.dir/data/raw
.dataset.list = "tcga_brca.RData"
#.dataset.list = "tcga_luad.RData"
#.dataset.list = "tcga_kirc.RData" 
#.dataset.list = "tcga_coad.RData" 

# choose a name for the current run
# all relative paths to project subfolders are set automatically
.current.run = "GeneExp_Met_2DS"
#.current.run = "LUAD_GeneExp_Met"
#.current.run = "KIRC_GeneExp_Met"
#.current.run = "COAD_GeneExp_Met"

source(file.path(.project.dir, "src/setDirPath.R"))

# Please edit the parameter files for synthetic and biological data
# The parameters are used for the current run

#####-------------------------------------------------------------#####
# install required packages                                           #
#####-------------------------------------------------------------#####

### Bioconductor
#source("https://bioconductor.org/biocLite.R")
#biocLite()
#biocLite("impute")                             # biological data, sCCA
#biocLite("org.Hs.eg.db")                       # comparison
#biocLite("GOstats")                            # comparison
#biocLite("graphite")                           # comparison
#biocLite("genefu")                             # comparison
#biocLite("SPIA")

### CRAN
#install.packages("httr")                       # TCGA Assembler
#install.packages("HGNChelper")                 # TCGA Assembler
#install.packages("RCurl")                      # TCGA Assembler
#install.packages("rjson")                      # TCGA Assembler
#install.packages("stringr")                    # TCGA Assembler
#install.pakcages("data.table")                 # TCGA Assembler 
#install.packages("gplots")                     # synthetic data
#install.packages("PMA")                        # sCCA
#install.packages("abind")                      # MALA
#install.packages("pROC")                       # drawROC
#install.packages("VennDiagram")                # comparison
#install.packages("xtable")                     # comparison
#install.packages("gridExtra")                  # comparison 
#install.packages("scales")                     # comparison
#install.packages("reshape2")                   # comparison
#install.packages("ggplot2")                    # comparison
#install.packages("Cairo")                      # Comparison

#####-------------------------------------------------------------#####
# run on biological datasets                                          #
#####-------------------------------------------------------------#####
# set parameteres as specified in the parameter file
source(file.path(.src.dir, "biologParameter.R"))

# starts TCGA data download using TCGA Assembler tool
source(file.path(.src.dir, "downloadTCGAData.R"))

# starts TCGA data preprocessing
source(file.path(.src.dir, "preprocTCGAData.R"))

# biological data exploration and transformation + sample reduction
source(file.path(.src.dir, "biologicalData.R"))

# do sCCA
source(file.path(.src.dir, "sCCA.R"))

# do NMF on biological data
source(file.path(.src.dir, "preprocForNMF.R"))
source(file.path(.src.dir, "NMF.R"))
source(file.path(.src.dir, "postprocOfNMF.R"))

# do pre- and postprocessing for MALA on biological datasets
# it is recommended to run MALA in a computationally more powerful
# linux environment due to the large size of the datasets
source(file.path(.src.dir, "preprocForMALA.R"))
#source(file.path(.src.dir, "MALA_linux.R")) # version to apply MALA 
                                             # to biologic data - run 
                                             # MALA only (not in project 
                                             # framework) on linux
source(file.path(.src.dir, "postprocOfMALA.R"))

# Compare result of each method (Venn diagrams, tables, ...)
source(file.path(.src.dir, "methodComparison.R"))

#####-------------------------------------------------------------#####
# run on synthetic datasets                                           #
#####-------------------------------------------------------------#####
# set parameteres as specified in the parameter file
source(file.path(.src.dir, "synthetParameter.R"))

# generate synthetic datasets 
source(file.path(.src.dir, "syntheticDataES.R"))

# do sCCA
source(file.path(.src.dir, "sCCA.R"))

# do NMF on synthetic data
source(file.path(.src.dir, "preprocForNMF.R"))
source(file.path(.src.dir, "NMF.R"))
source(file.path(.src.dir, "postprocOfNMF.R"))

# do pre- and postprocessing for MALA on synthetic datasets
# MALA is run in a computationally more powerful
# linux environment due to the large size of the datasets
source(file.path(.src.dir, "preprocForMALA.R"))
#source(file.path(.src.dir, "MALA.R")) # version to apply MALA to 
                                       # multiple synthetic datasets 
                                       # - run on a linux cluster
source(file.path(.src.dir, "postprocOfMALA.R"))

# Compare results of each method (Venn diagrams, boxplots, ...)
source(file.path(.src.dir, "synthetComparison.R"))








