#######################################################################
##### Project: IntMethodCompPublication                           #####
##### setDirPath.R                                                #####
#######################################################################

# main directory of analysis
.main.dir = .project.dir

# directory with binaries
.bin.dir = file.path(.main.dir, "bin")

# directory with source files
.src.dir = file.path(.main.dir, "src")

#directory with input data
.data.dir = file.path(.main.dir, "data")
.data.raw = file.path(.data.dir, "raw")

#directory with temporary data
.data.tmp = file.path(.data.dir, "tmp")

#directory with results
.results.dir = file.path(.main.dir, "results")

#----------------------------------------------------------------------

#directory with the results of the current run 
.run.dir = file.path(.results.dir, .current.run)
dir.create(path = .run.dir, recursive = T)


