multiNMF_mm = function(Xlist, K, maxiter, speak){
  #
  # Multiple NMF using euclidean distance update equations:
  #
  # Lee, D..D., and Seung, H.S., (2001), 'Algorithms for Non-negative Matrix
  # Factorization', Adv. Neural Info. Proc. Syst. 13, 556-562.
  #
  # INPUT: 
  # Xlist with i elements Xi
  # Xi (N,Mi): N (samples) x Mi (features) non negative input matrix
  # K        : Number of components
  # maxiter  : Maximum number of iterations to run
  # speak    : prints iteration count and changes in connectivity matrix 
  #            elements unless speak is 0
  #
  # OUTPUT:
  # Rlist containing W and i elements Hi:
  # W        : N x K matrix
  # Hi       : K x Mi matrix
  #
  # This script is based on TriNMF_mm.m function written by: 
  # Shihua Zhang
  # Computational Biology program in the Department of Biological Sciences
  # University of Southern California
  # zsh@amss.ac.cn
  # 2009/07/18
  
  ###############################
  # User adjustable parameters
  ###############################
  
  print_iter = 20 # iterations between print on screen and convergence test
  
  ############################################
  # test for negative values in input matrices
  ############################################
  
  if(min(unlist(lapply(Xlist, min))) < 0){
    cat("Input matrix elements can not be negative")
    return()
  }
  
  ################################################
  # test for same number of rows in input matrices
  ################################################
  
  ni = lapply(Xlist, nrow); mi = lapply(Xlist, ncol)

  if(length(unique(ni)) > 1){
    cat("Input matrices should have the same number of rows")
    return()
  }
  n = ni[[1]]
  
  #################################
  # initialize random W, H1 and H2
  #################################
  
  W = matrix(runif(n*K, min = 0, max = 1), nrow = n)
  Hlist = mapply(function(X,m) {matrix(runif(K*m, min=0, max=1), nrow=K, ncol=m)}, Xlist, mi)
  names(Hlist) = paste0("H", 1:length(Hlist))
  
  # use W*H to test for convergence
  Xr_old = lapply(Hlist, function(Hi,W) { W%*%Hi }, W)
 
  for(iter in 1:maxiter) {
    # Euclidean multiplicative method
    Hlist = mapply(function(Hi,Xi) { Hi*(t(W)%*%Xi)/((t(W)%*%W)%*%Hi + .Machine$double.eps) }, Hlist, Xlist) 
    w_num = W*t(Reduce("+", mapply(function(Hi,Xi) { Hi%*%Xi }, Hlist, lapply(Xlist, t), SIMPLIFY = F)))
    w_denom = W%*%(Reduce("+", lapply(Hlist, function(Hi) { Hi%*%t(Hi) }))) + .Machine$double.eps
    W = w_num/w_denom
    
    # print to screen
    if(iter %% print_iter == 0 & speak){
      Xr = mapply(function(Hi,W) { W%*%Hi }, Hlist, list(W))
      diff = Reduce("+", lapply(lapply(mapply("-", Xr_old, Xr), abs), sum))
      Xr_old = Xr
      eucl_dist = Reduce("+", mapply(nmf_euclidean_dist, Xlist, lapply(Hlist, function(Hi) { W%*%Hi })))
      errorx = Reduce("+", mapply(function(Xi, Hi) { mean(abs(Xi-W%*%Hi))/mean(Xi) }, Xlist, Hlist))
      cat("Iter = ", iter, " relative error = ", errorx, " diff = ", diff, 
          " eucl dist ", eucl_dist, "\n")
      
      if(errorx < 10^-5){
        break
      }
    }
  }
  return(c(W = list(W), Hlist))
}

nmf_euclidean_dist = function(X, Y){
  err = sum((X-Y)^2)
}

  
  


