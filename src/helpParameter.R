#######################################################################
##### file: helpParameter.R                                       #####
##### description: Extract the parameters used to create          #####
#####       synthetic datasets from their names.                  #####
##### input: dataset names                                        #####
##### author: B. Pucher                                           #####
#######################################################################

n_groups = sort(as.numeric(unlist(unique(strsplit(unlist(lapply(strsplit(datasets, split = "_"), "[", 3)), split = "gr")))))
p_de = sort(as.numeric(unlist(unique(strsplit(unlist(lapply(strsplit(datasets[grep("GeneExp", datasets)], split = "_"), "[", 4)), split = "p")))))
p_dm = sort(as.numeric(unlist(unique(strsplit(unlist(lapply(strsplit(datasets[grep("Methyl", datasets)], split = "_"), "[", 4)), split = "p")))))
c_de = sort(as.numeric(unlist(unique(strsplit(unlist(lapply(strsplit(datasets[grep("GeneExp", datasets)], split = "_"), "[", 5)), split = ".RData")))))
c_dm = sort(as.numeric(unlist(unique(strsplit(unlist(lapply(strsplit(datasets[grep("Methyl", datasets)], split = "_"), "[", 5)), split = ".RData")))))
  
