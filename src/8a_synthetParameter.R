#######################################################################
## Please define the parameters to be used for the current run       ##
## These are the parameters needed for the comparison of methods on  ##
##                                                                   ##
##  *** SYNTHETIC DATA ***                                           ##
##                                                                   ##
#######################################################################

# note: 'dataset' is data of one type e.g. gene expression data
#       'feature' is depending on data type a gene, a methyl site etc.

#####------------------------------------------------------------------
# parameters for generation of synthetic datasets
#####------------------------------------------------------------------

# number of samples in each condition of a dataset
# usually the sample groups should be of equal number
.n.tumor = 100
.n.normal = 100

# number of features in each simulated dataset type
.m.gene.exp = 1600
.m.methylation = 2400

# number of groups (pathways) preassumed in each dataset
.n.groups = 20

# proportion of pathways to be simulated as perturbed
.p.groups.d = 0.25

# percentage of features in a perturbed pathway to be simulated as different
# between conditions in each data type
# example: c(2, 10) means that two simulation runs are performed
#          one with 2 percent of features DE in perturbed pathways 
#          and another with 10 percent of features DE in perturbed pathways 
#          between conditions
.percentage.pw = c(10, 20, 50, 100)

# proportion of features additionally simulated as different between conditions
# and randomly spread on a dataset
# .p.random = 0.01
.p.random = 0

# Relative effect (peak) to add/subtract to the samples of one condition to simulate 
# DE features in a data type to induce weak, moderate or strong DE signals. 
# example: c(0.2, 0.4) means that two simulation runs are performed
#          one with a weak mean DE signal in DE features and 
#          and another with a moderate mean DE signal in DE features. 
.effect.size = c(0.2, 0.4, 0.8)

# number of datasets to create with each set of parameters (for performance statistics)
.stats.runs = 100
# .stats.runs = 10 # test the process first; the process works!

# create random seed for each run
set.seed(123456789)
.stats.seed = sample(111111111:888888888, .stats.runs, replace = FALSE)

# name of current synthetic data folder
gr = paste0(unique(c(min(.n.groups), max(.n.groups))), collapse = "_")
perc = paste0(unique(c(min(.percentage.pw), max(.percentage.pw))), collapse = "_")
es = paste0(unique(c(min(.effect.size), max(.effect.size))*100), collapse = "_")
.current.synthet = paste(
  "gr", gr,
  "perc", perc,
  "es", es,
  "bkgrnd", .p.random,
  sep = "_"
)

# specify data types to integrate
.data.to.integrate = c("GeneExp", "Methylation")

#####------------------------------------------------------------------
# parameters for sCCA on synthetic datasets
#####------------------------------------------------------------------
# iterations on optimization function
.niter = 50 # default 25

# increase e.g. when setting remain the same but synthetic data has changed
.serial = "ES"

# create directoy name part for current sCCA run
name = paste("exact", .serial, sep = "-")

# name of current sCCA result directory
.current.sCCA = paste("supervised", name, sep = "_")

#####------------------------------------------------------------------
# parameters for NMF on synthetic datasets
#####------------------------------------------------------------------

# number of basis vectors k to be used for NMF
# set this to the number of conditions in your data
.good.k = 2

# number of random initializations for the factor matrices W and Hi 
.nloop = 20

# number of iterations in each loop to minimize the reconstruction error
.maxiter = 1000

# increase e.g. when setting remain the same but synthetic data has changed
.serial = "ES"

# name of current NMF result directory
.current.NMF = paste0("supervised_k", .good.k, "-", .serial)    

#####------------------------------------------------------------------
# parameters for MALA on synthetic datasets
#####------------------------------------------------------------------

useMALA <- .useMALA_logi
if(useMALA){
  # number of cpus that can be dedicated to MALA
  .mala.cpus = 12
  
  # increase e.g. when setting remain the same but synthetic data has changed
  .serial = "ES"
  
  # parameters handed over to MALA
  # for documentation see the file bin\MALA\MALA-COMMAND-LINE-README.txt
  .sampling.type = 1
  .subsets = 1
  .set.covering.type = 2
  .beta = 50
  .graspsecs = 960 # default: 120, max: 960
  .graspiter = 100 # default: 100, max: 100000
  .numform = 3
  .exclusive.fs = 0
  
  # name of current MALA result folder
  .current.MALA = paste0(
    "ST", .sampling.type, "_SC", 
    .set.covering.type, "_B", .beta, "_Gs", 
    .graspsecs, "_Gi", .graspiter, "_NF", .numform,
    "_EF", .exclusive.fs, "-", .serial
  ) 
} else {
  .current.MALA <- NULL
}



#####------------------------------------------------------------------
# parameters for method comparison on synthetic data
#####------------------------------------------------------------------
.sCCA.run = .current.sCCA
.NMF.run = .current.NMF
.MALA.run = .current.MALA

add = ""
.current.comp = paste0(.current.synthet, add)

#####------------------------------------------------------------------
# do not change the .type parameter!
#####------------------------------------------------------------------
.type = "synthet"

#####------------------------------------------------------------------
# print out all hidden variables
#####------------------------------------------------------------------
ws = ls(all.names = TRUE)
hidden.vars = ws[grep("^\\.", ws)]
print(hidden.vars)

ls.str(all.names = TRUE)
