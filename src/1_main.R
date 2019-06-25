###############################################################################
##  Project: IntMethodCompPublication                                        ##
##  1_main.R                                                                 ##
##                                                                           ##
###############################################################################

rm(list = ls(all.names = TRUE))

# Please set path to the Project folder 
# .project.dir = "D:/Development/IntMethodCompPublication"
.project.dir = "C:/Users/gjo15/Documents/GitHub/IamComparison"

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

source(file.path(.project.dir, "src/2_setDirPath.R"))

# TCGA Assembler directory
# Please download this from https://github.com/compgenome365/TCGA-Assembler-2
.TCGA.assembler.dir = "F:/TCGA-Assembler-2/"

# Please edit the parameter files for synthetic and biological data
# The parameters are used for the current run

# If you execute this code on Windows or Macintosh machines, set this parameter
#   to FALSE
if(Sys.info()["sysname"] == "Linux"){
  .useMALA_logi <- TRUE
} else {
  .useMALA_logi <- FALSE
}

# 



######  install required packages  ############################################

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



######  run comparison on biological datasets  ################################
# set parameteres as specified in the parameter file
source(file.path(.src.dir, "3a_biologParameter.R"))

# starts TCGA data download using TCGA Assembler tool
source(file.path(.src.dir, "3b_downloadTCGAData.R"))

# starts TCGA data preprocessing
source(file.path(.src.dir, "3c_preprocTCGAData.R"))

# biological data exploration and transformation + sample reduction
source(file.path(.src.dir, "3d_biologicalData.R"))

# do sCCA
source(file.path(.src.dir, "4_sCCA.R"))

# do NMF on biological data
source(file.path(.src.dir, "5a_preprocForNMF.R"))
source(file.path(.src.dir, "5b_NMF.R"))
source(file.path(.src.dir, "5c_postprocOfNMF.R"))

# do pre- and postprocessing for MALA on biological datasets
# it is recommended to run MALA in a computationally more powerful
# linux environment due to the large size of the datasets
source(file.path(.src.dir, "6a_preprocForMALA.R"))
# source(file.path(.src.dir, "6b_MALA_linux.R"))
# # version to apply MALA to biologic data - run MALA only (not in project 
# #   framework) on linux

source(file.path(.src.dir, "6c_postprocOfMALA.R"))

# Compare result of each method (Venn diagrams, tables, ...)
source(file.path(.src.dir, "7_methodComparison.R"))



######  run comparison on synthetic datasets  #################################

# set parameteres as specified in the parameter file
source(file.path(.src.dir, "8a_synthetParameter.R"))

# generate synthetic datasets 
system.time(
  source(file.path(.src.dir, "8b_syntheticDataES.R"))
)
# about 30-ish minutes(?) for 10 repititions. This script cleans the global
#   environment before execution, so "a" was removed. I'm going to delete these
#   results.
# 24.03667, 25.44217 min for 10 reps
# The gene expression data is 100 x 1600; the methlyation data is 100 x 2400
# 53.41383 min for 100 replicates
# For 4 x 3 design, this takes 67.2555 min for 100 reps
# For 16 design points and 100 replicates: 86.22083 min



# do sCCA
system.time(
  source(file.path(.src.dir, "4_sCCA.R"))
)
# 17.8418 hrs for 10 repititions (9 design points). 36.19844 hrs for 15 reps by
#   12 design points.
# 140.3611 hours for 43 reps by 16 design points

# do NMF on synthetic data
system.time(
  source(file.path(.src.dir, "5a_preprocForNMF.R"))
)
# 7.067167 min (9 design points). 8.718 min for 15 reps by 12 design points.
# 31.7375, 30.79133 min for 43 reps by 16 design points

system.time(
  source(file.path(.src.dir, "5b_NMF.R"))
)
# 15.59804 hrs for 10 repititions (9 design points). 33.98197 hrs for 15 reps
#   by 12 design points
# 123.162 hours for 43 reps by 16 design points

system.time(
  source(file.path(.src.dir, "5c_postprocOfNMF.R"))
)
# 3.747 min (9 design points).


useMALA <- .useMALA_logi
if(useMALA){
  # do pre- and postprocessing for MALA on synthetic datasets
  # MALA is run in a computationally more powerful
  # linux environment due to the large size of the datasets
  source(file.path(.src.dir, "6a_preprocForMALA.R"))
  # 9.077 min
  
  source(file.path(.src.dir, "6b_MALA.R"))
  # version to apply MALA to multiple synthetic datasets - run on a linux
  #   cluster
  
  source(file.path(.src.dir, "6c_postprocOfMALA.R"))
}


# Compare results of each method (Venn diagrams, boxplots, ...)
source(file.path(.src.dir, "9_synthetComparison.R"))




