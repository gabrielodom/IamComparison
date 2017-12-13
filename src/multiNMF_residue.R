multiNMF_residue = function(Xlist, kstart, kend, nloop, verbose) {
  #
  # Model Selection for NMF ()
  #
  # Xlist: Xi (N,Mi): N (samples) x Mi (features) non negative input matrix
  #
  # kstart, kend : range of values of k to test consensus for, in order 
  #                to find a good k
  #
  # nloop: number of initial conditions per k
  #        (start with 10 or 20, check with more)
  # 
  # verbose : prints iteration count and changes in connectivity matrix 
  #           elements if not set to 0
  #
  #
  # consensus : 3d array of consensus matrices 
  #             dimensions : kend x M x M 
  #             Values in consensus(1:kstart-1,:,:) should be ignored
  #
  # This script is based on TriNMF_residue.m function written by: 
  # Shihua Zhang
  # Computational Biology program in the Department of Biological Sciences
  # University of Southern California
  # zsh@amss.ac.cn
  # 2009/07/18
  #
  
  # test for negative values in v 
  if(min(unlist(lapply(Xlist, min))) < 0) {
    cat("matrix entries can not be negative")
    return()
  } 
  
  if(min(unlist(lapply(Xlist, rowSums))) == 0) {
    cat("not all entries in a row can be zero")
    return()
  }
  
  Obj_residue = lapply(Xlist, function(x) numeric(kend))
  
  for(j in kstart:kend) {
    if (verbose){
      cat("rank ", j, "\n")
    }
    
    objhistory = lapply(Xlist, function(x) numeric(nloop))
    
    for(iloop in 1:nloop){
      if(verbose){
        cat("iteration ", iloop, "\n")
      }
      
      maxiter = 200; speak = 1
      WHlist = multiNMF_mm(Xlist, j, maxiter, speak)
      
      # compute residue
      newobj = mapply(function(X, H) {sum((X-WHlist$W%*%H)^2)}, Xlist, WHlist[grep("H",names(WHlist))])
      
      objhistory = mapply(function(oh, no) {oh[iloop] = no; oh}, objhistory, newobj, SIMPLIFY = FALSE)
      
    }
    
    Obj_residue = mapply(function(or, oh) {or[j] = oh; or}, 
                         Obj_residue, lapply(objhistory, min), SIMPLIFY = F)
    
  }
  
  return(Obj_residue)
}
