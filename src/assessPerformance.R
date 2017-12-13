library("gplots")

#####------------------------------------------------------------------
# Assess classification power of methods
#####------------------------------------------------------------------
assessPerformance = function(gene.set, method, biolog.data, dist.fun, 
                             hclust.fun, merge = FALSE, plot = FALSE){
  # subset of data with features selected by method
  data.subset = mapply(function(d, r) {
    d[which(rownames(d) %in% r),]
  }, d = biolog.data, r = gene.set, SIMPLIFY = FALSE)
  
  # create RowSideColors vector 
  ft.type = cbind(c(rep(names(data.subset), unlist(lapply(data.subset, nrow)))), "merged")
  row.vc = as.character(factor(ft.type[,1], labels = rainbow(length(unique(ft.type[,1])), 
                                                             start = 0.3, end = 0.7)))
  
  if(merge == TRUE){
    data.subset = list(do.call(rbind, data.subset))
    names(data.subset) = "merged"
  }
  
  # mask clustering functions with different methods
  dist.mask = function(x) dist(x, method = dist.fun)
  hclust.mask = function(x) hclust(x, method = hclust.fun) 
  
  
  # visualize classification in heatmap and save to file 
  cfmtx = invisible(mapply(function(x, n) {
    cl = hclust.mask(dist.mask(t(x)))
    tr = cutree(cl, k = 2)
    tr.type = sub("_.*", "", colnames(x), fixed = FALSE)
    col.vc = as.character(factor(tr.type, labels = c("red", "black")))

    if(plot == TRUE){
      # draw and save heatmap
      row.vc = row.vc[which(ft.type == n, arr.ind = TRUE)[,"row"]]  
      png(file.path(sub.dir.figures, paste(method, n, hclust.fun, "HM.png", sep = "_")),
      width = 900, height = 900)
      hm = heatmap.2(x, distfun = dist.mask, hclustfun = hclust.mask,
                     trace = "none", scale = "row", col = redblue(12),
                     #breaks = seq(-6, 6, length.out = 13),
                     ColSideColors = col.vc, RowSideColors = row.vc,
                     labRow = "", labCol = "",
                     margins = c(1,1), dendrogram = "column",
                     main = paste0("Sample Clustering by ", method, 
                                   " (", n, ", ", hclust.fun, ")"),
                     keysize = 1)
      # title(paste0("Sample Clustering by ", method, " (", n, ", ", hclust.fun, ")"), 
      #       cex.main = 2, outer = TRUE)
      legend("topleft", unique(tr.type), fill = unique(col.vc), 
             inset = c(0, 0.21), title = "true class", cex = 1)
      legend("topleft", unique(ft.type[which(ft.type == n, arr.ind = TRUE)[,"row"],1]), 
             fill = unique(row.vc), inset = c(0, 0.35), title = "data type", cex = 1)

      dev.off()
    }    

    
    # create and print confusion matrix
    class = c("a", "b")
    cfmtx = matrix(NA, ncol = 2, nrow = 2)
    dimnames(cfmtx) = list(paste(method, class, sep = "."), 
                           c("biolog.tu", "biolog.no"))
    #cfmtx = table(tr, tr.type)
    #rownames(cfmtx) = paste(method, unique(tr.type), sep = "_")
    
    for(i in c("tumor", "normal")){
      for(j in c(1,2)){
        cfmtx[grep(class[j], rownames(cfmtx)), grep(substr(i, 1,2), colnames(cfmtx))] =
          length(which(tr.type == i & tr == j))
      }
    }
    if(cfmtx[1,1] > cfmtx[1,2] & cfmtx[2,1] < cfmtx[2,2]){
      rownames(cfmtx) = sub("a", "tu", rownames(cfmtx))
      rownames(cfmtx) = sub("b", "no", rownames(cfmtx))
    } else {
      rownames(cfmtx) = sub("a", "no", rownames(cfmtx))
      rownames(cfmtx) = sub("b", "tu", rownames(cfmtx))
    }
    
    return(list(cfmtx))
    
  }, x = data.subset, n = names(data.subset)))  
  
  return(cfmtx)
}

getAccuracy = function(cfmtx){

  TP = cfmtx[grep("tu", rownames(cfmtx)), grep("tu", colnames(cfmtx))]
  FP = cfmtx[grep("tu", rownames(cfmtx)), grep("no", colnames(cfmtx))]
  TN = cfmtx[grep("no", rownames(cfmtx)), grep("no", colnames(cfmtx))]
  FN = cfmtx[grep("no", rownames(cfmtx)), grep("tu", colnames(cfmtx))]
  
  acc = (TP+TN)/(TP+FN+TN+FP)
  
  return(acc)
  
}

getTPR = function(cfmtx){
  TP = cfmtx[grep("tu", rownames(cfmtx)), grep("tu", colnames(cfmtx))]
  FN = cfmtx[grep("no", rownames(cfmtx)), grep("tu", colnames(cfmtx))]
  TPR = TP/(TP+FN)
  return(TPR)
}

getFPR = function(cfmtx){
  FP = cfmtx[grep("tu", rownames(cfmtx)), grep("no", colnames(cfmtx))]
  TN = cfmtx[grep("no", rownames(cfmtx)), grep("no", colnames(cfmtx))]
  FPR = FP/(FP+TN)
  return(FPR)
}
