# Job 1: NNMF
# Gabriel Odom
# 2019-07-02

# I found out about jobs in RStudio. Apparently, I can have a script run in the
#   background of an R session, instead of having two windows for the same
#   project open.


######  Setup  ################################################################
rm(list = ls(all.names = TRUE))

# Please set path to the Project folder 
.project.dir <- "C:/Users/gjo15/Documents/GitHub/IamComparison"
.dataset.list <- "tcga_brca.RData"
.current.run <- "GeneExp_Met_2DS"

source(file.path(.project.dir, "src/2_setDirPath.R"))

# TCGA Assembler directory
# Please download this from https://github.com/compgenome365/TCGA-Assembler-2
.TCGA.assembler.dir = "F:/TCGA-Assembler-2/"

if(Sys.info()["sysname"] == "Linux"){
  .useMALA_logi <- TRUE
} else {
  .useMALA_logi <- FALSE
}



######  run comparison on synthetic datasets  #################################

# set parameteres as specified in the parameter file
source(file.path(.src.dir, "8a_synthetParameter.R"))

# NMF
source(file.path(.src.dir, "5a_preprocForNMF.R"))
source(file.path(.src.dir, "5b_NMF.R"))
source(file.path(.src.dir, "5c_postprocOfNMF.R"))
