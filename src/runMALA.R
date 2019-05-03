#######################################################################
##### file: runMALA.R                                             #####
##### description: This file is supposed to be executed with      #####
#####       RScript in file MALA.R to run MALA on Linux.          #####
##### input: in-/output paths as arguments according to MALA.R    #####
##### output: output files of MALA in corresponding out directory #####                          
##### author: B. Pucher                                           #####
##### date created: 18/05/2017                                    #####
##### last change:  09/11/2017                                    #####
#######################################################################
args = commandArgs(trailingOnly = TRUE)

des = as.numeric(args[1])
mala.tmp.sub = args[2]
file.name = args[3]
mala.res.dir = args[4]
max.splits = as.numeric(args[5])
run.seed = as.numeric(args[6])

src.dir = sub(paste0(.Platform$file.sep, "results", .Platform$file.sep, ".*"), 
              paste0(.Platform$file.sep, "src", .Platform$filesep), mala.res.dir)

source(file.path(src.dir, "extractMALAFeatures.R"))
invisible(source(file.path(src.dir, "8a_synthetParameter.R")))

setwd(mala.tmp.sub)

log.con = file("runMALALog.txt")
sink(file = log.con, type = "output")
flush(log.con)

invisible(file.create("timing.txt"))

cat(paste("-----", Sys.time(), "-----"), "\n")  
cat("Desired number of features:", des, "\n")

res.dir = rev(unlist(strsplit(mala.res.dir, split = .Platform$file.sep)))[1]
sub.dir.files =  sub(res.dir, "", mala.res.dir)
features = character(0); s = 1

set.seed(run.seed)
while(length(features) < des & s <= max.splits){
  # permute features in input file and save to TMP-folder

  if(.set.covering.type == 1){
    input.data = read.csv(file.name, header = FALSE, stringsAsFactors = FALSE)
    input.data.perm = input.data[c(1:2, sample(3:nrow(input.data), nrow(input.data)-2, 
                                               replace = FALSE)),]
    write.table(input.data.perm, file = file.name, sep  = ",",
                quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  cat("Set MALA parameters. \n")
  cat("MALA seed: ", run.seed, "\n")  
  
  # overwrite MALA parameters file for windows and linux
  parameters = file.path(mala.tmp.sub, "parameters.dat")
  
  writeLines(c(
    paste("PADDING", 0, sep = "\t"), 
    paste("FASTAINPUT",	0, sep = "\t"),
    paste("CSVINPUT",	1, sep = "\t"),
    paste("CLUSTERING", 0, sep = "\t"),
    paste("SAMPLINGTYPE",	.sampling.type, sep = "\t"),
    paste("SUBSETS",	.subsets, sep = "\t"),
    paste("PERCSLICING",	80, sep = "\t"),
    paste("MAXSPLITS",	8, sep = "\t"),
    paste("ENTROPY",	0.6, sep = "\t"),
    paste("SCTYPE",	.set.covering.type, sep = "\t"),
    paste("BETA", 	.beta, sep = "\t"),
    paste("GRASPSECS",	.graspsecs, sep = "\t"),
    paste("GRASPITER",	.graspiter, sep = "\t"),
    paste("GRASPALPHA",	0.6, sep = "\t"),
    paste("RANDOMSEED",	run.seed, sep = "\t"),
    paste("POSCOST",	1, sep = "\t"),
    paste("NEGCOST",	1, sep = "\t"),
    paste("NUMFORM",	.numform, sep = "\t"),
    paste("EXCLUSIVEFS", .exclusive.fs, sep = "\t")),
    con = parameters)

  if(file.exists("outputMicroA")){
    # delete existing output folder in MALA tmp dir
    status_del = file.remove(dir("outputMicroA", full.names = TRUE))
  }
  
  cat("Start MALA on dataset", file.name, "\n")
  cmd.mala.linux = paste("/usr/bin/time -v ./runMala.sh", file.name, "2>> timing.txt")
  system(cmd.mala.linux, wait = TRUE)

  # rename result files after current run and copy them to results folder 
  invisible(file.rename(dir("outputMicroA", full.names = TRUE),
                        sub("DMB_s0", paste0("DMB_s", s-1), 
                            dir("outputMicroA", full.names = TRUE))))
  
  # remove temporary results because they are VERY big
  system("rm -r TMP_*")
  
  # copy result files to MALA results dir of current dataset 
  status_cp_res = file.copy(from = dir("outputMicroA", full.names = T),
                            to = mala.res.dir, overwrite = TRUE)
  
  # extract formulas and compare number of features to desired
  features = unique(unlist(extractMALAFeatures(res.dir, desired = .m.gene.exp+.m.methylation)))
  cat("Number of features:", length(features), "\n")
  
  s = s+1
  run.seed = sample(111111111:888888888, 1)
}

cat(paste("-----", Sys.time(), "-----"), "\n")
sink()
close(log.con)

# copy parameter file to output folder
status_cp_param = file.copy(from = "parameters.dat", 
                            to = file.path(mala.res.dir, "parameters.dat"))

# copy log file to output folder
status_cp_log = file.copy(from = "runMALALog.txt", 
                          to = file.path(mala.res.dir, "runMALALog.txt"))

# copy time-log file to output folder
status_cp_time = file.copy(from = "timing.txt",
                          to = file.path(mala.res.dir, "timing.txt"))






