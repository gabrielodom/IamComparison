source(file.path(.src.dir, "assessPerformance.R"))

#####------------------------------------------------------------------
# Split datasets k times randomly into training and test set of samples. 
#####------------------------------------------------------------------
randomSplit = function(data, subsets, p.train, 
                         samples.paired = FALSE, fname.add = ""){

  cat("Start preparing dataset", fname.add ,"for Cross validation ...\n")
  
  for(k in 1:subsets){

    # Split datasets 
    tu.idx = sample(grep("tu", colnames(data[[1]])), 
                    round(length(grep("tu", colnames(data[[1]])))*p.train/100))
    if(samples.paired){
      # select the normal sample originating from the same patient as the tumor sample
      cat("Assuming paired tumor and normal tissue samples. \n")
      no.idx = grep("no", colnames(data[[1]]))[tu.idx]
    } else {
      no.idx = sample(grep("no", colnames(data[[1]])),
                      round(length(grep("no", colnames(data[[1]])))*p.train/100))
    }
    
    train.idx = c(tu.idx, no.idx)
    test.idx = setdiff(c(1:ncol(data[[1]])), train.idx)
    data.train = lapply(data, function(x) x[,train.idx]); names(data.train) = names(data)
    data.test = lapply(data, function(x) x[,test.idx]); names(data.test) = names(data)
    saveRDS(data.train, file = file.path(sub.dir.RData, paste0("data_train", k, fname.add, ".RData")))
    saveRDS(data.test, file = file.path(sub.dir.RData, paste0("data_test", k, fname.add, ".RData")))
  }
  
  cat("Dataset", k, "times randomly split into training and test data.\n")
}


#####------------------------------------------------------------------
# Calculate confusion matrix and accuracy of classification power of 
# method on training and test dataset 
#####------------------------------------------------------------------
evaluatePerformance = function(method, subsets, fname.add = ""){
  
  acc.train.vc = numeric(subsets)
  acc.test.vc = numeric(subsets)
  cf.train.list = list()
  cf.test.list = list()
  TPR.train.vc = numeric(subsets)
  TPR.test.vc = numeric(subsets)
  FPR.train.vc = numeric(subsets)
  FPR.test.vc = numeric(subsets)

  for(k in 1:subsets){
    
    data.train = readRDS(file.path(.run.dir, data.dir, "RData", 
                                   paste0("data_train", k, fname.add, ".RData")))
    data.test = readRDS(file.path(.run.dir, data.dir, "RData",
                                  paste0("data_test", k, fname.add, ".RData")))
    result.train = readRDS(file.path(sub.dir.RData, paste0(method, "_result_train", k, fname.add, ".RData")))
    
    # calculate confusion matrix for current cross validation run
    cfmtx.train = assessPerformance(result.train, method = paste0(method, "_tr", k),
                                    data.train, dist.fun = "euclidean", 
                                    hclust.fun = "ward.D2", merge = TRUE, plot = TRUE)
    cfmtx.test = assessPerformance(result.train, method = paste0(method, "_te", k), 
                                   data.test, dist.fun = "euclidean",
                                   hclust.fun = "ward.D2", merge = TRUE, plot = TRUE)
    
    acc.train = lapply(cfmtx.train, getAccuracy)
    acc.test = lapply(cfmtx.test, getAccuracy)
    
    acc.train.vc[k] = mean(do.call(c, acc.train)); names(acc.train.vc)[k] = paste0("train",k)
    acc.test.vc[k] = mean(do.call(c, acc.test)); names(acc.test.vc)[k] = paste0("test",k)
    
    TPR.train = lapply(cfmtx.train, getTPR)
    TPR.test = lapply(cfmtx.test, getTPR)
    
    TPR.train.vc[k] = mean(do.call(c, TPR.train)); names(TPR.train.vc)[k] = paste0("train",k)
    TPR.test.vc[k] = mean(do.call(c, TPR.test)); names(TPR.test.vc)[k] = paste0("test",k)
    
    FPR.train = lapply(cfmtx.train, getFPR)
    FPR.test = lapply(cfmtx.test, getFPR)
    
    FPR.train.vc[k] = mean(do.call(c, FPR.train)); names(FPR.train.vc)[k] = paste0("train",k)
    FPR.test.vc[k] = mean(do.call(c, FPR.test)); names(FPR.test.vc)[k] = paste0("test",k)
    
    cf.train.list[[k]] = cfmtx.train; names(cf.train.list)[k] = paste0("train", k)
    cf.test.list[[k]] = cfmtx.test; names(cf.test.list)[k] = paste0("test", k)
    
    cat("Confusion matrix for training set", k, "\n")
    print(cfmtx.train)
    cat("Classification accuracy on training set", k, ":", 
        paste(names(acc.train), unlist(lapply(acc.train, round, 3)) * 100, "%"), "\n\n") 
    
    cat("Confusion matrix of test set: \n")
    print(cfmtx.test)
    cat("Classification accuracy on test set:", 
        paste(names(acc.test), unlist(lapply(acc.test, round, 3)) * 100, "%"), "\n\n") 
  }  
  cat("Average accuracy of", method, "on training datasets ( splits =",
      .subsets.val, "):", round(mean(acc.train.vc), 3) * 100, "%  (SD:", 
      round(sd(acc.train.vc), 3) * 100,  "%) \n")
  cat("Average accuracy of", method, "on test datasets ( splits =",
      .subsets.val, "):", round(mean(acc.test.vc), 3) * 100, "%  (SD:",
      round(sd(acc.test.vc), 3) * 100, "%) \n")
  
  cat("Average TPR of", method, "on training datasets (splits =",
      .subsets.val, "):", round(mean(TPR.train.vc), 3) * 100, "%  (SD:",
      round(sd(TPR.train.vc), 3) * 100, "%) \n")
  cat("Average TPR of", method, "on test datasets (splits =",
      .subsets.val, "):", round(mean(TPR.test.vc), 3) * 100, "%  (SD:",
      round(sd(TPR.test.vc), 3) * 100, "%) \n")
  
  cat("Average FPR of", method, "on training datasets (splits =",
      .subsets.val, "):", round(mean(FPR.train.vc), 3) * 100, "%  (SD:",
      round(sd(FPR.train.vc), 3) * 100, "%) \n")
  cat("Average FPR of", method, "on test datasets (splits =",
      .subsets.val, "):", round(mean(FPR.test.vc), 3) * 100, "%  (SD:",
      round(sd(FPR.test.vc), 3) * 100, "%) \n")

  saveRDS(acc.train.vc, file = file.path(sub.dir.RData, paste0(method, "_acc_train", fname.add, ".RData")))
  saveRDS(acc.test.vc, file = file.path(sub.dir.RData, paste0(method, "_acc_test", fname.add, ".RData")))
  saveRDS(TPR.train.vc, file = file.path(sub.dir.RData, paste0(method, "_TPR_train", fname.add, ".RData")))
  saveRDS(TPR.test.vc, file = file.path(sub.dir.RData, paste0(method, "_TPR_test", fname.add, ".RData")))
  saveRDS(FPR.train.vc, file = file.path(sub.dir.RData, paste0(method, "_FPR_train", fname.add, ".RData")))
  saveRDS(FPR.test.vc, file = file.path(sub.dir.RData, paste0(method, "_FPR_test", fname.add, ".RData")))
  saveRDS(cf.train.list, file = file.path(sub.dir.RData, paste0(method, "_cfm_train", fname.add, ".RData")))
  saveRDS(cf.test.list, file = file.path(sub.dir.RData, paste0(method, "_cfm_test", fname.add, ".RData")))

}


