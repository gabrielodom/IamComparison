#######################################################################
## Please define the parameters to be used for the current run       ##
## These are the parameters needed for the comparison of methods on  ##
##                                                                   ##
##  *** BIOLOGICAL DATA ***                                          ##
##                                                                   ##
#######################################################################

# note: 'dataset' is data of one type e.g. gene expression data
#       'feature' is depending on data type a gene, a methylon site etc.

#####------------------------------------------------------------------
# parameters required by TCGA Assembler 
#####------------------------------------------------------------------
.TCGA.assembler.dir = "D:/TCGA-Assembler/"
.sCancer = "BRCA"

#####------------------------------------------------------------------
# parameters for the preprocessing of biological dataset
#####------------------------------------------------------------------

# specify data types to integrate
.data.to.integrate = c("GeneExp", "Methylation")

# name of directory to save preprocessed data
.current.biolog = "less_10p_NA-eval-seed2"

# Features to select in percent. (Refere to log file in biologicalData!)
.pr = c(1.14, 1.36) # BRCA
#.pr = c(1.95, 0.56) # LUAD
#.pr = c(3.53, 1.9) # KIRC
#.pr = c(1.57, 3.22) # COAD
  
#####------------------------------------------------------------------
# Number of subsets and percentage of training samples for validation
#####------------------------------------------------------------------
# should data be prepared for cross validation? 
.prepare.CV = TRUE

.subsets.val = 10
.pr.train = 80

# should cross validatin be performed or skipped? 
.do.CV = TRUE

# seed for random percentage split of samples
.biolog.seed = 123456789

# seed for method runs on cross validation subsets
set.seed(987654321)
.biolog.seed.cv = sample(111111111:888888888, .subsets.val, replace = FALSE) 
  
#####------------------------------------------------------------------
# parameters for sCCA on biological datasets
#####------------------------------------------------------------------
# iterations on optimization function
.niter = 100 # default 25

# create directoy name part for current sCCA run
name = paste0(sub("\\.", "_", sum(.pr)))

# set e.g. when settings remain the same but data set has changed
serial = ""

# name of current sCCA result directory
.current.sCCA = paste0(.current.biolog, "-", name, serial)

#####------------------------------------------------------------------
# parameters for MALA on biological datasets
#####------------------------------------------------------------------

# set e.g. when settings remain the same but data set has changed
serial = ""

# parameters handed over to MALA
# for documentation see the file bin\MALA\MALA-COMMAND-LINE-README.txt
.sampling.type = 1 # 1 = random percentage split, 2 = cross validation
.subsets = 100
.set.covering.type = 2 # 1 = linear, 2 = quadratic
.beta = 50
.graspsecs = 960
.graspiter = 100
.numform = 3
.exclusive.fs = 0

# name of current MALA result folder
.current.MALA = paste0("ST", .sampling.type, "_SC", 
                       .set.covering.type, "_B", .beta, "_Gs", 
                       .graspsecs, "_Gi", .graspiter, "_NF", .numform,
                       "_EF", .exclusive.fs, "-", .current.biolog, serial)

#####------------------------------------------------------------------
# parameters for NMF on biological datasets
#####------------------------------------------------------------------

# seed for NMF on biological data
.biolog.NMF.seed = .biolog.seed

# number of basis vectors k to be used for NMF
# set this to the number of conditions in your data
.good.k = 2

# number of random initializations for the factor matrices W and Hi  
.nloop = 50

# number of iterations in each loop to minimize the reconstruction error
.maxiter = 1000

# set e.g. when settings remain the same but data set has changed
serial = ""

# name of current NMF result directory
.current.NMF = paste0(.current.biolog, "-", sub(".", "_", sum(.pr), fixed = T), serial)

#####------------------------------------------------------------------
# parameters for method comparison on biological data
#####------------------------------------------------------------------
.sCCA.run = .current.sCCA
.NMF.run = .current.NMF
.MALA.run = .current.MALA

add = "-Gs960"
.current.comp = paste0(.current.biolog, "-", 
                       sub(".", "_", sum(.pr), fixed = T), add)

#####------------------------------------------------------------------
# do not change the .type parameter!
#####------------------------------------------------------------------
.type = "biolog"

#####------------------------------------------------------------------
# print out all hidden variables
#####------------------------------------------------------------------
ws = ls(all.names = TRUE)
hidden.vars = ws[grep("^\\.", ws)]
print(hidden.vars)

