multiNMF_residue_comodule = function(Xlist, K, nloop, maxiter){
  
  verbose = 1
  # nloop = 100
  # maxiter = 1000
  
  ni = lapply(Xlist, nrow); mi = lapply(Xlist, ncol)
  n = ni[[1]]
  
  bestW = matrix(0, nrow = n, ncol = K)
  bestHlist = mapply(function(X,m) { matrix(0, nrow=K, ncol=m)}, Xlist, mi)
  names(bestHlist) = paste0("H", 1:length(bestHlist))
  
  bestobj = rep(1000000000, length(bestHlist))
  
  for (iloop in 1:nloop) {
    if(verbose) {
      cat("iteration ", iloop, "\n")
    }
    
    speak = 1
    WHlist = multiNMF_mm(Xlist, K, maxiter, speak)
    
    # compute residue
    newobj = mapply(function(X,H) { sum((X-WHlist$W%*%H)^2) },
                    Xlist, WHlist[grep("H", names(WHlist))])
    
    if (any(mapply("<", newobj, bestobj))) {
      bestobj = newobj
      bestW = WHlist$W
      bestHlist = WHlist[grep("H", names(WHlist))]
      
    }
  }
  
  W = bestW; Hlist = bestHlist
  
  return(c(W = list(W), Hlist))
}


