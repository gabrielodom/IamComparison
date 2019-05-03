# Run MALA subsequently on all .csv files in this folder.
# Input files must be named <input_filename>_DMB.csv
# Results are saved in folders named <input_filename>_out

# MALA bin files, MALA_linux.R and the input files in csv format 
# must be located in the same folder and this folder must be the 
# R working directory (navigate to the folder and then start R).

log.con = file("MALALog.txt")
sink(file = log.con, type = "output")
flush(log.con)


inputs = list.files(pattern = "\\.csv$")
if(length(inputs) == 0){
  cat("No input files found. Please provide csv files to analyze in",
      "directory:", getwd(), "\n")
} else {

invisible(file.create("timing.txt"))
  
  for(i in 1:length(inputs)){
    
    # delete files in existing output folder outputMicroA
    cat("Deleting output files of previous run of MALA ... \n")
    status_del = file.remove(dir("outputMicroA", full.names = T))
    print(table(status_del))
    
    cat("\n")
    cat("Current file: ", inputs[i], "\n\n")
    
    # run MALA
    cmd = paste0("/usr/bin/time -v ./runMala.sh ", inputs[i], " 2>> timing.txt")
    cat("MALA started ... \n")
    system(cmd)
    
    # copy parameters file to results folder
    cat("MALA finished. \nCopy parameters file to output folder ... \n")
    status_cp_param = file.copy(from = "parameters.dat",
                                to = file.path("outputMicroA", "parameters.dat"))
    cat(status_cp_param, "\n")
    
  	# copy time-log file to output folder
  	cat("Copy timing-log file to output folder ... \n")
	  status_cp_time = file.copy(from = "timing.txt",
                          to = file.path("outputMicroA", "timing.txt"))
						  
    # rename result folder
    cat("Create new folder and move result files there ... \n")
    res.dir = sub("DMB.csv", "out", inputs[i])
    dir.create(res.dir)
    status_cp_res = file.copy(from = dir("outputMicroA", full.names = T),
                              to = res.dir, overwrite = T)
    print(table(status_cp_res))
    
  }
}


sink()
close(log.con)
