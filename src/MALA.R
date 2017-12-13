#######################################################################
##### file: MALA.R                                                #####
##### input: synthetic datasets created in a number of runs       #####
##### output: logic formulas for sample classification .csv       #####
##### packages: --                                                #####
##### author: B. Pucher                                           #####  
##### date created: 12/04/2017                                    #####
##### last change:  09/11/2017                                    #####
#######################################################################
rm(list = ls())

#######################################################################
#####                                                             #####
##### FUNCTIONS                                                   #####
#####                                                             #####
#######################################################################


#######################################################################
#####                                                             #####
##### MAIN SECTION                                                #####
#####                                                             #####
#######################################################################
current.run = .current.MALA
synthet.run = .current.synthet
sub.dir.name = file.path(paste0(.type, "MALA"), current.run)
source(file.path(.src.dir, "setSubDirPath.R"))

log.con = file(file.path(sub.dir.files, "MALALog.txt"))
sink(file = log.con, type = "output")
flush(log.con)

cat(paste("-----", Sys.time(), "-----"), "\n")

mala.bin.dir = file.path(.bin.dir, "MALA-UBUNTU-DEBIAN-LINUX-synthet")

# create directory for results of current MALA run 
mala.tmp.dir = file.path(.data.tmp, paste0("TMP-", current.run))
dir.create(mala.tmp.dir, showWarnings = FALSE)

if(.type == "synthet"){

  # function definitions
  source(file.path(.src.dir, "extractMALAFeatures.R"))
  
  stats.runs = .stats.runs
  cpus = .mala.cpus; cat("Number of screens used for MALA:", cpus, "\n")
  sub.dir.files.tmp = sub.dir.files
  run = 1

  # get dataset names in one stats run
  datasets = list.files(file.path(.run.dir, "syntheticData", synthet.run, "RData",
                                  paste0("run", run)))
  datasets = datasets[-grep("_Diff|test|train|PWs", datasets)]

  # get n_groups, p_de, p_dm, c_de, c_dm from datasets names
  source(file.path(.src.dir, "helpParameter.R"))
  
  for(run in 1:stats.runs){
    
    sub.dir.files = file.path(sub.dir.files.tmp, paste0("run", run))
    if(length(list.files(sub.dir.files, pattern = "_out")) > 0){
      # delete existing output in MALA results dir
      cat("Delete results of previous executions of MALA on the same dataset.\n")
      system(paste("rm -r", file.path(sub.dir.files, "*_out")))
    }

    cat("Create temporary directory for MALA stats run number ", run, ".\n")
    mala.stats.dir = file.path(mala.tmp.dir, paste0("run", run))
    dir.create(mala.stats.dir, showWarnings = FALSE)

    run.seed = .stats.seed[run]
    cat("Seed used:", run.seed, "\n")
          
    for(k in 1:length(n_groups)){
      for(i in 1:length(p_de)){
        for(j in 1:length(c_de)){
          
          fname_ge = paste("GeneExp_synth_", n_groups[k], "gr_", p_de[i], "p_", c_de[j], sep = "")
          fname_met = paste("Methyl_synth_",  n_groups[k], "gr_", p_dm[i], "p_", c_dm[j], sep = "")
          params = paste0(letters[k], letters[i], letters[j])
          file.name = paste(.type, params, "DMB.csv", sep = "_")
          mala.input = file.path(sub.dir.files, file.name)
          
          cat("Create temporary directory for MALA run on dataset: ", file.name, "\n")
          mala.tmp.sub = file.path(mala.stats.dir, paste0("TMP-", sub("\\.csv", "", file.name)))
          dir.create(mala.tmp.sub, showWarnings = FALSE)
          cat("Results will be temporarily saved in", mala.tmp.sub, "\n")
          
          # copy bin files required by MALA and current dataset to temp dir
          cat("Copy files to temporary directory. \n")
          invisible(file.copy(from = dir(mala.bin.dir, full.names = T),
                              to = mala.tmp.sub, overwrite = TRUE))
          invisible(file.copy(from = mala.input, to = mala.tmp.sub,
                              overwrite = TRUE))
          
          # directory where MALA results of current dataset shoud be saved to
          mala.res.dir = file.path(sub.dir.files, paste(.type, params, "out", sep = "_"))
          dir.create(mala.res.dir, showWarnings = FALSE)
          
          # desired number of features to select
          desired = sum(unlist(getDesiredSynthet(.current.synthet, run)))
          
          # start MALA runs on current synthetic dataset
          max.splits = 100
          
          # allow a maximum number of cpu screens but one screen per cpu at most 
          # for parallel execution of MALA on linux cluster
          wait = TRUE
          while(wait == TRUE){ 
            
            # count screens
            sockets = system("screen -ls | grep Socket", intern = TRUE)
            screens = as.numeric(unlist(strsplit(sockets, split = " "))[1])
            
            if(screens < min(30, cpus)){

              # Start new instance of MALA 
              cat("Start MALA in new screen ... \n")
              system(paste("screen -md -S", paste0(params, run), "Rscript", 
                           file.path(.src.dir, "runMALA.R"), desired, 
                           mala.tmp.sub, file.name, mala.res.dir, max.splits, run.seed))

              # wait a few seconds to shift/delay next start of MALA
              Sys.sleep(10)

              cat("Results are saved in ", mala.res.dir, "\nMALA has finished!\n")
              
              wait = FALSE
 
            } else {
              Sys.sleep(10)
            }
          }
          
        }
      }
    }
  }
  
} else {
  print("Cannot execute MALA on biological dataset.\n")
}


cat(paste("-----", Sys.time(), "-----"), "\n")
cat("Parameters used and results of MALA saved to results folder. \n")

sink()
close(log.con)
    
 






