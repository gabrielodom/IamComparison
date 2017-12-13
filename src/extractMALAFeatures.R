#######################################################################
##### file: extractMALAFeatures.R                                 #####
##### description: Extract feature names from logic formulas      #####
#####       resulting from MALA                                   #####
##### input: path to folder containing MALA formula files         #####
##### output: list of features selected in each dataset           #####                          
##### author: B. Pucher                                           #####
##### date created: 11/05/2017                                    #####
##### last change:  09/11/2017                                    #####
#######################################################################

#####------------------------------------------------------------------
# Extract feature names from logic formulas resulting from MALA.
#####------------------------------------------------------------------
extractMALAFeatures = function(mala.out.dir, desired){
  
  # read in all lines in the formula files resulting from MALA
  mala.raw = lapply(file.path(sub.dir.files, mala.out.dir, 
                              list.files(file.path(sub.dir.files, mala.out.dir), "formulas.csv")),
                    readLines)
  
  cat("Grep result files containing formulas and extract gene list... \n")
  cat("Current folder:", mala.out.dir, "\n")
  
  # extract lines containing the actual logic formulas
  # code adapted from Chris Fischer - start
  to.remove = "^CL.+|^OR|^Cov.+|^ Fal.+|^ Sco.+|^ FP.+"
  mala.formulas = lapply(mala.raw, function(x) x[grep(to.remove, x, invert = TRUE)])
  mala.formulas = lapply(mala.formulas, function(x) x[which(x != "")])
  
  feature.names = lapply(mala.formulas, function(x) { strsplit(
    #remove commas and parentheses
    gsub("\\,|\\(|\\)","",
         #remove numbers and relational operators ahead of gene names
         gsub("[0-9\\.-]+[<>=]+","",
              #remove points and numbers following a relational operator 
              gsub("[<>=]+[0-9\\.-]+","", x))),
    #split literals and clauses to obtain gene symbols
    split=" +AND +| +OR +") })
  
  gene.names.ge.flat = character(0)
  gene.names.ge = lapply(feature.names, function(x) {
    tmp = unique(unlist(
      lapply(mapply("[", x,
                    lapply(x, grep, pattern = "^ge")),
             gsub, pattern = "ge.", replacement = "", fixed = T)))
    # get empty character instead of NULL if no features are detected
    if(length(tmp) == 0){ return("") } else { return(tmp) }
  })
  
  gene.names.met.flat = character(0)
  gene.names.met = lapply(feature.names, function(x) {
    tmp = unique(unlist(
      lapply(mapply("[", x,
                    lapply(x, grep, pattern = "^met")),
             gsub, pattern = "met.", replacement = "", fixed = T)))
    # get empty character instead of NULL if no features are detected
    if(length(tmp) == 0){ return("") } else { return(tmp) }
  })

  # code adapted from Chris Fischer - end
  
  gene.names.flat = character(0)
  i = 1; while(length(gene.names.flat) <  desired & i <= length(mala.raw)){
    gene.names.ge.flat = union(gene.names.ge.flat, gene.names.ge[[i]])
    gene.names.met.flat = union(gene.names.met.flat, gene.names.met[[i]])
    gene.names.flat = c(gene.names.ge.flat, gene.names.met.flat)
    i = i+1
  }
  
  cat("Results of ", i-1, " MALA run(s) accumulated to get close to the desired",
      "number of ", desired, " features.\n")
  
  genes.list = list(gene.names.ge.flat, gene.names.met.flat)
  names(genes.list) = c("GeneExp", "Methyl")
  
  return(genes.list)
}

#####------------------------------------------------------------------
# Get desired nb of features to be selected by MALA in synthetic data
#####------------------------------------------------------------------
getDesiredSynthet = function(synthet.run = "", run = ""){
  
  # POSITIVES: differentially expressed/methylated features
  pos_ge = readRDS(file.path(.run.dir, "syntheticData", synthet.run, 
                             "RData", paste0("run", run), 
                             paste0(fname_ge, "_Diff.RData")))
  pos_met = readRDS(file.path(.run.dir, "syntheticData", synthet.run, 
                              "RData", paste0("run", run),
                              paste0(fname_met, "_Diff.RData")))
  
  desired = list(length(pos_ge), length(pos_met))
  
  return(desired)
}

#####------------------------------------------------------------------
# Get desired nb of features to be selected by MALA in biological data
#####------------------------------------------------------------------
getDesiredBiolog = function(data.raw, data.names){
  
  # Number of assumed POSITIVES:
  desired = mapply(function(d, p) {
    round(dim(d)[1] * (p/100))
  }, d = data.raw, p = .pr, SIMPLIFY = FALSE)
  
  names(desired) = data.names
  
  return(desired)
}
