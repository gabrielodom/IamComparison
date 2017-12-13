#######################################################################
##### Project: IntMethodCompPublication                           #####
##### setSubDirPath.R                                             #####
#######################################################################

#directory with the results of current sub procedure
sub.dir = file.path(.run.dir, sub.dir.name)

#if(!file.exists(sub.dir)){
  dir.create(path = sub.dir, recursive = T)
  
  #directory with the file results of the current run
  sub.dir.files = file.path(sub.dir, "files")
  dir.create(sub.dir.files, recursive = T)
  
  #directory with the RData results of the current run
  sub.dir.RData = file.path(sub.dir,"RData")
  dir.create(sub.dir.RData, recursive = T)
  
  #directory with the figure results of the current run
  sub.dir.figures = file.path(sub.dir,"figures")
  dir.create(sub.dir.figures, recursive = T)
#}


