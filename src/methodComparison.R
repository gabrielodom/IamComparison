#######################################################################
##### file: methodComparison.R                                    #####
##### input: flat gene lists resulting from methods to compare    #####
##### output: Venn diagrams, tables, p-values of overlaps and/or  #####
##### individual results of GO and pathway analysis.              #####
##### packages: GOstats, org.Hs.eg.db, graphite, genefu, gplots   ##### 
#####           VennDiagram, xtable, gridExtra                    #####
##### author: B. Pucher                                           #####  
##### date created: 30/04/2015                                    #####
##### last change:  24/01/2018                                    #####
#######################################################################

rm(list = ls())

library(org.Hs.eg.db)
library(GOstats)
library(graphite)
library(genefu)
library(SPIA)

# CRAN
library(VennDiagram)
library(xtable)
library(gridExtra)
library(gplots)
library(ggplot2)
library(reshape2)
library(scales)

#######################################################################
#####                                                             #####
##### FUNCTIONS                                                   #####   
#####                                                             #####
#######################################################################

#####------------------------------------------------------------------
# get corresponding Entrez IDs for given gene symbols
#####------------------------------------------------------------------
mapSymbol2EG = function(symbols, get.from.raw.data = TRUE, named = FALSE){

  if(get.from.raw.data) { # pick EntrezIDs from input dataset
 
    load(file.path(.data.raw, .dataset.list)) # inputDataList
    sy.EG = inputDataList[[1]]$Des
    if(named == FALSE) {
      entrez.genes = unname(unlist(lapply(symbols, function(s) {
        sy.EG[which(sy.EG[,"GeneSymbol"] == s), "EntrezID"]
      })))
    } else {
      entrez.genes = unlist(lapply(symbols, function(s) {
        id = sy.EG[which(sy.EG[,"GeneSymbol"] == s), "EntrezID"]
        id = ifelse(length(id) == 0, yes = NA, no = id)
        names(id) = s; id
        }))
    }

    return(entrez.genes)
    
  }  else {

    y = org.Hs.egSYMBOL2EG
    mapped.genes = mappedkeys(y) # Ids mapped to a gene symbol
    yy = as.list(y[mapped.genes])
    
    zz = as.list(org.Hs.egALIAS2EG) # Ids mapped to several symbols
    zz = zz[!is.na(zz)]
    
    entrez.genes = yy[symbols]
    nn.symbols = setdiff(symbols, names(entrez.genes)[which(!is.na(
      names(entrez.genes)))])
    alias.EG = zz[nn.symbols]
    
    # using unlist, a sequential number is added to the name of list 
    # elements containing more than one element !!!!
    entrez.genes.compl = unique(unlist(c(entrez.genes, alias.EG)[which(!is.na(names(c(
      entrez.genes, alias.EG))))]))
    
    return(entrez.genes.compl) 
    
  }
}

#####------------------------------------------------------------------
# get corresponding gene names for given Entrez IDs
#####------------------------------------------------------------------
mapEG2Name = function(EG){

  x = org.Hs.egGENENAME
  mapped.genes = mappedkeys(x)
  
  xx = as.list(x[mapped.genes])
  names = xx[EG]
  
  return(names)

}

#####------------------------------------------------------------------
# get Entrez ID for feautre identifier using lookup table
#####------------------------------------------------------------------
mapFeaturesToEG = function(features, des, ft.col) {
  
  ft.row = which(des[,ft.col] %in% features)
  ft.EG = des[ft.row, "EntrezID"]
  
  return(ft.EG)
}

#####------------------------------------------------------------------
# get gene symbol for feautre identifier using lookup table
#####------------------------------------------------------------------
mapFeaturesToGenes = function(features, des, ft.col) {
  
  ft.row = which(des[,ft.col] %in% features)
  ft.genes = des[ft.row, "GeneSymbol"]
  
  return(ft.genes)
}

#####------------------------------------------------------------------
# get gene symbols missing in lookup table from Entrez ID
#####------------------------------------------------------------------
mapFeaturesToGenesHs = function(features.EG){
  
  mapping = select(org.Hs.eg.db, keys = features.EG, columns = "SYMBOL",
                   keytype = "ENTREZID")
  
  return(mapping$SYMBOL)
}

#####------------------------------------------------------------------
# create parameter objects for GO analysis 
#####------------------------------------------------------------------
createParams = function(geneset, universe, ontology, 
                        param.conditional = FALSE, param.pval = 0.05){

  # Remove Genes which are not contained in the universe
  geneset.un = geneset[which(is.element(geneset, universe))]

  params = new("GOHyperGParams", ontology = ontology,
               geneIds = geneset.un, universeGeneIds = universe,
               annotation = "org.Hs.eg.db", pvalueCutoff = param.pval, 
               testDirection = "over", conditional = param.conditional)
  
  return(params)
}

#####------------------------------------------------------------------
# Get genes in cancer signature with name "cancer.set".
#####------------------------------------------------------------------
getCancerGenes = function(cancer.set) {
  
  cancer.genes.EG = 
    switch(cancer.set,
           kegg = {
             pic = pathways("hsapiens", "kegg")$'Pathways in cancer'
             unique(c(pic@edges$src, pic@edges$dest))
           },
           pam = {
             # require(genefu)
             as.character(unique(na.omit(pam50$centroids.map$EntrezGene.ID)))
           },
           sorlie = {
             # require(genefu)
             as.character(unique(na.omit(ssp2003$centroids.map$EntrezGene.ID)))
           },
           hu = {
             # require(genefu)
             as.character(unique(na.omit(ssp2006$centroids.map$EntrezGene.ID)))
           },
           bedo = {
             temp = c("PDCD1", "TBX21", "IDOX", "IFNG", "CCL5", 
                      "GZMB", "PRF1", "GNLY", "CD274", "CTLA4", 
                      "CD19", "CXCL9", "CD8B", "STAT1", "CXCL13",
                      "IRF1", "CXCL10", "IGKC", "IL12B", "FOXP3", 
                      "IL12A")
             unique(mapSymbol2EG(temp, get.from.raw.data = FALSE))
           },
           immune = {
             temp = read.table(file.path(.data.dir, "databases", "immune_genes.txt"),
                               header = TRUE, stringsAsFactors = FALSE)
             as.character(unique(temp$EntrezGene))
           },
           target = {
             temp = read.table(file.path(.data.dir, "databases", "target_db_annotation.txt"),
                               header = TRUE, stringsAsFactors = FALSE)
             as.character(unique(temp$EntrezID))
           })
  
  return(cancer.genes.EG)
  
}

#####------------------------------------------------------------------
# Intersection of selected genes with known cancer genesets
#####------------------------------------------------------------------
intersectCancerGenes = function(symbols, cancer.set){
  
  cancer.genes.EG = getCancerGenes(cancer.set)
  genes.EG = mapSymbol2EG(symbols)
  intersection = intersect(genes.EG, cancer.genes.EG)
  
  return(intersection) 
}

#####------------------------------------------------------------------
# Assess significance of overlap btw. method results on feature level
# and save histogram 
#####------------------------------------------------------------------
assessSignOverlap = function(ft.set.list, gene.set.list, ft.universe, 
                             nsample = 1000){
  
  overlap = Reduce(intersect, gene.set.list)
  
  noverlap = numeric(nsample)
  for(i in 1:nsample) {
    sample1 = mapply(function(f, u) { sample(size = length(f), x = u)},
                     f = ft.set.list[[1]], u = ft.universe )
    sample2 = mapply(function(f, u) { sample(size = length(f), x = u)},
                     f = ft.set.list[[2]], u = ft.universe )
    sample3 = mapply(function(f, u) { sample(size = length(f), x = u)},
                     f = ft.set.list[[3]], u = ft.universe )
    
    sample1.total = unique(unlist(sample1))
    sample2.total = unique(unlist(sample2))
    sample3.total = unique(unlist(sample3))
    
    noverlap[i] = length(Reduce(intersect, list(sample1.total, 
                                                sample2.total,
                                                sample3.total)))
  }
  
  hist(noverlap, prob = T, 
       xlab = "Size of Overlap", xlim = c(0, length(overlap)),
       breaks = seq(min(noverlap), max(noverlap), by = 1),
       main = paste0("Distribution of Overlap btw. Method Results"))
  axis(side = 1, at = c(length(overlap)))
  lines(density(noverlap, adjust = 2), lwd = 2, col = "darkblue")
  abline(v = length(overlap))
  text(x = length(overlap), y = max(density(noverlap, adjust = 2)$y)*0.8, 
       labels = paste("overlap"), pos = 2)
  invisible(dev.print(png, file.path(sub.dir.figures, paste0("Signf_overlap_hist.png")),
                      width = 1000))
  
  Fn = ecdf(noverlap)  
  p.greater = 1-Fn(length(overlap))
  return(round(p.greater, 2))
}

#####------------------------------------------------------------------
# Assess significance of overlap btw. method results and cancer 
# genes and save histograms
#####------------------------------------------------------------------
assessSign = function(gene.set, cancer.set, name, overlap, universe, nsample = 20000) {
  
  cancer.genes.EG = getCancerGenes(cancer.set)
  
  in.uni = intersect(cancer.genes.EG, data.all.universe.EG)
    
  noverlap = numeric(nsample)
  for(i in 1:nsample) {
    gs.sample = sample(size = length(gene.set), x = universe)
    noverlap[i] = length(intersect(gs.sample, in.uni))
  }
  
  hist(noverlap, prob = T, breaks = seq(-0.5, max(noverlap)+0.5, by = 1), 
       xlab = "Size of Overlap [count]", 
       main = paste0(strsplit(name, split = "_")[[1]][2], ": Distribution of Overlap with ", 
                     toupper(cancer.set), " Cancer Genes"))
  axis(side = 1, at = c(length(overlap)))
  lines(density(noverlap, adjust = 5), lwd = 2, col = "darkblue")
  abline(v = length(overlap))
  text(x = length(overlap), y = max(density(noverlap, adjust = 2)$y)*0.9, 
       labels = paste(strsplit(name, split = "_")[[1]][2], "overlap"), pos = 4)
  invisible(dev.print(png, file.path(sub.dir.figures, paste0(name, "_overlap_hist.png")),
                      width = 1000))
  
  Fn = ecdf(noverlap)  
  p.greater = 1-Fn(length(overlap))
  return(round(p.greater, 2))
  
}

#####------------------------------------------------------------------
# Test for over-representation of cancer genes in method
# results with Fisher's exact test. 
#####------------------------------------------------------------------
assessOdds = function(gene.set.EG, cancer.set, universe){
  
  cancer.genes.EG = getCancerGenes(cancer.set)
  
  a = intersect(cancer.genes.EG, gene.set.EG)
  b = setdiff(gene.set.EG, cancer.genes.EG)
  c = intersect(setdiff(universe, gene.set.EG), cancer.genes.EG)
  d = intersect(setdiff(universe, cancer.genes.EG), setdiff(universe, gene.set.EG))
  
  ct = matrix(c(length(a), length(b), length(c), length(d)),
              nrow = 2, byrow = FALSE,
              dimnames = list(Signature = c("Yes", "No"),
                              Method = c("Yes", "No")))
  
  fp = fisher.test(ct, alternative = "greater")$p.value
  
  return(round(fp,5))
}

#####------------------------------------------------------------------
# Create tables of a given data frame in txt, tex and pdf format
#####------------------------------------------------------------------
createTables = function(df, fname){

  # table in text format
  write.table(df, file = file.path(sub.dir.files, paste0(fname, ".txt")),
              sep = "\t", quote = F, row.names = F)
  
  # table in LaTeX format
  xdf = xtable(df)
  print(xdf, file = file.path(sub.dir.files, paste0(fname, ".tex")),
        sanitize.text.function=function(x){x},
        include.rownames = F)
  
  # table as picture 
  pdf(file.path(sub.dir.files, paste0(fname, ".pdf")), 
      width = 12, height = max(11, nrow(df)/2.5-4))
  grid.newpage()
  grid.table(df, rows = NULL)
  invisible(dev.off())
  
}

#####------------------------------------------------------------------
# Merge tables of features, GO terms or RP resulting from each method.
#####------------------------------------------------------------------
mergeTables = function(fname, header, merge_by){

  # create list with tables to merge 
  methods_tab = lapply(methods, function(m) {
    if(file.exists(file.path(sub.dir.files, paste0(m, "_", fname)))){
      read.delim(file.path(sub.dir.files, paste0(m, "_", fname)), 
                 header = TRUE, as.is = TRUE) 
    } else {
      data.frame(matrix(vector(), 0, length(header), dimnames = list(c(), header)),
                 stringsAsFactors = FALSE)
    } })
  names(methods_tab) = methods
  
  # merge tables and write to file
  merged_tab = Reduce(function(x, y) merge(x, y, by = merge_by, sort = FALSE, all = TRUE), 
                      methods_tab)
  colnames(merged_tab)[which(!is.element(colnames(merged_tab), merge_by))] = 
          rep(methods, each = (ncol(merged_tab)-length(merge_by))/length(methods))
  write.table(merged_tab[order(merged_tab[,grep("Size$|GeneSymbol", colnames(merged_tab))], 
                               decreasing = TRUE),], file.path(sub.dir.files, fname),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
}

#####------------------------------------------------------------------
# Calculate log Fold Change for given expression profiles of two groups
# with paired samples
#####------------------------------------------------------------------
getLFC = function(data, groups, use.diff = FALSE){
  
  if(use.diff == TRUE){
    
    case = data[,grep(groups[1], colnames(data))]
    control = data[,grep(groups[2], colnames(data))]
    diff = rowMeans(case - control)
    names(diff) = rownames(data)
    return(diff)
    
  } else {
    
    case = 2^data[,grep(groups[1], colnames(data))]
    control = 2^data[,grep(groups[2], colnames(data))]
    lfc = rowMeans(log2(case/control))
    names(lfc) = rownames(data)
    return(lfc)
  }
}

#####------------------------------------------------------------------
# Load performance measures from methods result directory. 
#####------------------------------------------------------------------
loadPerfMeasure = function(pm){
  
  sCCA_meas_train = readRDS(file.path(.run.dir, "biologsCCA", .sCCA.run, 
                                      "RData", paste0("sCCA_", pm, "_train.RData")))
  sCCA_meas_test = readRDS(file.path(.run.dir, "biologsCCA", .sCCA.run, 
                                     "RData", paste0("sCCA_", pm, "_test.RData")))
  NMF_meas_train = readRDS(file.path(.run.dir, "biologNMF", .NMF.run, 
                                     "RData", paste0("NMF_", pm, "_train.RData")))
  NMF_meas_test = readRDS(file.path(.run.dir, "biologNMF", .NMF.run, 
                                    "RData", paste0("NMF_", pm, "_test.RData")))
  MALA_meas_train = readRDS(file.path(.run.dir, "biologMALA", .MALA.run, 
                                      "RData", paste0("MALA_", pm, "_train.RData")))
  MALA_meas_test = readRDS(file.path(.run.dir, "biologMALA", .MALA.run, 
                                     "RData", paste0("MALA_", pm, "_test.RData")))
  
  meas_arr = do.call(rbind, list(sCCA_tr = sCCA_meas_train, sCCA_te = sCCA_meas_test,
                                 NMF_tr = NMF_meas_train, NMF_te = NMF_meas_test,
                                 MALA_tr = MALA_meas_train, MALA_te = MALA_meas_test))
  
  return(meas_arr) 
}

#####------------------------------------------------------------------
# Create Boxplot of performance measures and save to file.
#####------------------------------------------------------------------
plotPerfMeasureBP = function(perf_ggplot, fname = "Boxplot_performance.eps"){
  
  pr = ggplot(data = perf_ggplot, aes(x = measure, y = value, 
                                      fill = Method, col = Method)) +
    geom_boxplot(position = position_dodge(0.85)) +
    scale_fill_manual(values = brewer_pal(palette = "Pastel2")(3)) +
    scale_color_manual(values = brewer_pal(palette = "Set2")(3)) +
    scale_x_discrete(labels = levels(perf_ggplot$measure)) +
    scale_y_continuous(breaks = seq(0, 110, 20),
                       limits = c(0,110), oob = rescale_none) +
    facet_grid(~set, scales = "free") +
    labs(list(x = "Performance measure", y = "%")) +
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          panel.grid.major.y = element_line(colour = "grey50", linetype = 2),
          panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(),
          text = element_text(size = 25), axis.text = element_text(size = 25),
          legend.position = "top", legend.background = element_rect(colour = "white"),
          legend.direction = "horizontal", legend.key = element_rect(fill = "white"),
          legend.text = element_text(size = 20))
  
  postscript(file = file.path(sub.dir.figures, fname), 
             onefile = FALSE, width = 12.5, height = 5.5, paper = "special", 
             family = "serif", horizontal = FALSE)
  print(pr)
  dev.off()
  
}

#######################################################################
#####                                                             #####
##### MAIN SECTION                                                #####
#####                                                             #####
#######################################################################
current.run = .current.comp
sub.dir.name = file.path("biologComparison", current.run)
source(file.path(.src.dir, "setSubDirPath.R"))
source(file.path(.src.dir, "assessPerformance.R"))

set.seed(.biolog.seed)

log.con = file(file.path(sub.dir.files, "methodComparisonLog.txt"))
sink(file = log.con, type = "output")
flush(log.con)

cat("Compare results of method runs:\n sCCA:", .sCCA.run, "\n NMF:", 
    .NMF.run, "\n MALA:", .MALA.run, "\n")

#####------------------------------------------------------------------
# load lists of features resulting from the methods to compare
# results of sCCA, NMF and MALA
#####------------------------------------------------------------------
sCCA_result = readRDS(file.path(.run.dir, "biologsCCA", .sCCA.run, 
                                "RData", "sCCA_result.RData"))
load(file.path(.run.dir, "biologNMF", .NMF.run, "RData", 
               "NMF_result.RData"))
MALA_result = readRDS(file.path(.run.dir, "biologMALA", .MALA.run, "RData", 
                                "MALA_result.RData"))

# Numbers of unique features resulting from each method per data type 
cat("Size of feature sets resulting from sCCA run:", 
    mapply(paste, names(sCCA_result), lapply(sCCA_result, length)), "\n")
cat("Size of feature sets resulting from NMF run:", 
    mapply(paste, names(NMF_result), lapply(NMF_result, length)), "\n")
cat("Size of feature sets resulting from MALA run:", 
    mapply(paste, names(MALA_result), lapply(MALA_result, length)), "\n")

#####------------------------------------------------------------------
# map feature lists to genes (symbols, Entrez IDs) to compare
# results of sCCA, NMF and MALA
#####------------------------------------------------------------------
# Get gene symbols from feature identifiers (use lookup table)
feature_des = lapply(file.path(.run.dir, "preprocTCGAData-2", "RData", 
                                paste(.data.to.integrate, "Des.RData", sep = "_")), 
                     readRDS)
names(feature_des) = .data.to.integrate

# Try to find missing Gene symbols by Entrez ID
feature_des$GeneExp[which(feature_des$GeneExp[,"GeneSymbol"] == "?"), "GeneSymbol"] <-
  mapFeaturesToGenesHs(feature_des$GeneExp[which(feature_des$GeneExp[,"GeneSymbol"] == "?"), "EntrezID"])

sCCA_result_genes_maps = mapply(mapFeaturesToGenes, features = sCCA_result,
                           des = feature_des, ft.col = c("EntrezID", "REF"))
NMF_result_genes_maps = mapply(mapFeaturesToGenes, features = NMF_result,
                          des = feature_des, ft.col = c("EntrezID", "REF"))
MALA_result_genes_maps = mapply(mapFeaturesToGenes, features = MALA_result,
                           des = feature_des, ft.col = c("EntrezID", "REF"))   

cat("\n")

# Numbers of maps to genes (symbols) in each method per data type
cat("Number of maps to genes resulting from sCCA run:", 
    mapply(paste, names(sCCA_result_genes_maps), 
           lapply(sCCA_result_genes_maps, length)), "\n")
cat("Number of maps to genes resulting from NMF run:", 
    mapply(paste, names(NMF_result_genes_maps), 
           lapply(NMF_result_genes_maps, length)), "\n")
cat("Number of maps to genes resulting from MALA run:", 
    mapply(paste, names(MALA_result_genes_maps), 
           lapply(MALA_result_genes_maps, length)), "\n")

cat("\n")

# Get unique gene symbols resulting from methods
sCCA_result_genes = lapply(sCCA_result_genes_maps, 
                           function(x) unique(x[which(!is.na(x))]))
NMF_result_genes = lapply(NMF_result_genes_maps, 
                          function(x) unique(x[which(!is.na(x))]))
MALA_result_genes = lapply(MALA_result_genes_maps, 
                           function(x) unique(x[which(!is.na(x))]))

# Numbers of unique genes (symbols) resulting from each method per data type
cat("Number of unique genes (symbols) resulting from sCCA run:", 
    mapply(paste, names(sCCA_result_genes), 
           lapply(sCCA_result_genes, length)), "\n")
cat("Number of unique genes (symbols) resulting from NMF run:", 
    mapply(paste, names(NMF_result_genes), 
           lapply(NMF_result_genes, length)), "\n")
cat("Number of unique genes (symbols) resulting from MALA run:", 
    mapply(paste, names(MALA_result_genes), 
           lapply(MALA_result_genes, length)), "\n")

rm(sCCA_result_genes, NMF_result_genes, MALA_result_genes)

cat("\n")

# Numbers of features which can not be mapped to a gene symbol
cat("Number of features resulting from sCCA lacking a gene symbol:",
    mapply(paste, names(sCCA_result_genes_maps),
           lapply(sCCA_result_genes_maps, 
                  function(x) length(which(is.na(x))))), "\n")
cat("Number of features resulting from NMF lacking a gene symbol:",
    mapply(paste, names(NMF_result_genes_maps),
           lapply(NMF_result_genes_maps, 
                  function(x) length(which(is.na(x))))), "\n")
cat("Number of features resulting from MALA lacking a gene symbol:",
    mapply(paste, names(MALA_result_genes_maps),
           lapply(MALA_result_genes_maps, 
                  function(x) length(which(is.na(x))))), "\n")

cat("\n")

# Add Entrez ID column to Methylation lookup table 
# Map symbols to Entrez IDs from GeneExp datasets:
methyl.EG = mapSymbol2EG(feature_des$Methylation[,"GeneSymbol"], named = TRUE)
methyl.EG.na = names(methyl.EG[which(is.na(methyl.EG))])
# Map symbols not found in GeneExp dataset
methyl.na.unique = unique(methyl.EG.na[which(!is.na(methyl.EG.na))])
methyl.EG.na.Hs = select(org.Hs.eg.db, keys = methyl.na.unique,
                      columns = "ENTREZID", keytype = "SYMBOL")
# Merge symbol mappings
for(sym in methyl.na.unique) {
  id = methyl.EG.na.Hs[which(methyl.EG.na.Hs$SYMBOL == sym), "ENTREZID"]
  methyl.EG[sym] = ifelse(length(id) == 0, yes = NA, no = id)
}
# Add symbol mappings to annotation table
methyl.EG.len = ifelse(names(methyl.EG)[length(methyl.EG)] == "", 
                       yes = length(methyl.EG)-1, no = length(methyl.EG)) 
# identical(names(methyl.EG)[1:methyl.EG.len], feature_des$Methylation[,"GeneSymbol"]) # TRUE
feature_des$Methylation = cbind(feature_des$Methylation, EntrezID = methyl.EG[1:methyl.EG.len])

# Numbers of maps to genes (Entrez ID) in each method per data type (use lookup table)
sCCA_result_genes_maps_EG = mapply(mapFeaturesToEG, features = sCCA_result,
                              des = feature_des, ft.col = c("EntrezID", "REF"))
NMF_result_genes_maps_EG = mapply(mapFeaturesToEG, features = NMF_result,
                             des = feature_des, ft.col = c("EntrezID", "REF"))
MALA_result_genes_maps_EG = mapply(mapFeaturesToEG, features = MALA_result,
                              des = feature_des, ft.col = c("EntrezID", "REF"))   

# Get unique genes (Entrez IDs) resulting from methods
sCCA_result_genes_EG = lapply(sCCA_result_genes_maps_EG, 
                           function(x) unique(x[which(!is.na(x))]))
NMF_result_genes_EG = lapply(NMF_result_genes_maps_EG, 
                          function(x) unique(x[which(!is.na(x))]))
MALA_result_genes_EG = lapply(MALA_result_genes_maps_EG, 
                           function(x) unique(x[which(!is.na(x))]))

# Numbers of genes (Entrez IDs) resulting from each method per data type
cat("Number of unique genes (Entrez IDs) resulting from sCCA run:", 
    mapply(paste, names(sCCA_result_genes_EG), 
           lapply(sCCA_result_genes_EG, length)), "\n")
cat("Number of unique genes (Entrez IDs) resulting from NMF run:", 
    mapply(paste, names(NMF_result_genes_EG), 
           lapply(NMF_result_genes_EG, length)), "\n")
cat("Number of unique genes (Entrez IDs) resulting from MALA run:", 
    mapply(paste, names(MALA_result_genes_EG), 
           lapply(MALA_result_genes_EG, length)), "\n")
cat("\n")

methods = c("sCCA", "NMF", "MALA")

#####------------------------------------------------------------------
# comparison on gene level
#####------------------------------------------------------------------

cat("Compare results on feature level:\n")
cat("Plotting Venn-Diagrams ...\n")

# overlap of features resulting from each dataset 
for(i in 1:length(.data.to.integrate)) {
  sets.list = list(sCCA_result[[i]], NMF_result[[i]], MALA_result[[i]])
  names(sets.list) = paste0(methods, " (", lapply(sets.list, function(x) {
    length(unique(x)) }), ")")
  tmp = venn.diagram(sets.list, 
               filename = file.path(sub.dir.figures, paste0("venn_features_", 
                  .data.to.integrate[i], ".png")),
               main = paste("Feature sets:", .data.to.integrate[i]), 
               main.cex = 2, cex = 1.6,
               cat.cex = 2, cat.dist = c(0.07, 0.07, 0.05), 
               print.mode = c("raw"), sigdigs = 1,
               lty = "solid", cat.pos = c(-30,30,180), 
               euler.d = FALSE, scale = FALSE)
  
}

cat("Compare results on gene level:\n")
cat("Plotting Venn-Diagrams ...\n")

# overlap of genes resulting from each dataset 
for(i in 1:length(.data.to.integrate)) {
  sets.list = list(sCCA_result_genes_EG[[i]], NMF_result_genes_EG[[i]], MALA_result_genes_EG[[i]])
  names(sets.list) = paste0(methods, " (", lapply(sets.list, length), ")")
  tmp = venn.diagram(sets.list, 
               filename = file.path(sub.dir.figures, 
                                    paste0("venn_genes_", .data.to.integrate[i], ".png")),
               main = paste("Gene sets:", .data.to.integrate[i]), 
               main.cex = 2, cex = 1.6,
               cat.cex = 2, cat.dist = c(0.07, 0.07, 0.05), 
               print.mode = c("raw"), sigdigs = 1,
               lty = "solid", cat.pos = c(-30,30,180), 
               euler.d = FALSE, scale = FALSE)
  
}

# overlap of features (accumulated)
sCCA_result_total = unique(unlist(sCCA_result))
NMF_result_total = unique(unlist(NMF_result))
MALA_result_total = unique(unlist(MALA_result))

ft.set.list = list(sCCA_result, NMF_result, MALA_result)
names(ft.set.list) = methods
ft.set.list.total = list(sCCA_result_total, NMF_result_total, MALA_result_total)
names(ft.set.list.total) = paste0(methods, " (", lapply(ft.set.list.total, length), ")")
 
cat("Length of resulting feature sets:", names(ft.set.list.total), "\n")
cat("Cardinality of the union of feature sets:", 
    length(Reduce(union, ft.set.list.total)), "\n")
invisible(venn.diagram(ft.set.list.total, 
             filename = file.path(sub.dir.figures, "venn_features_all.png"),
             main = "Feature sets: all data types", main.cex = 2,
             cex = 1.6, cat.cex = 2, cat.pos = c(-25, 25, 180),
             cat.dist = c(0.06, 0.06, 0.04),
             print.mode = c("raw"), sigdigs = 1,
             lty = "solid", euler.d = FALSE, scale = FALSE))

tmp <- venn.diagram(ft.set.list.total, width = 7000, height = 7000,
                    units = "px", resolution = 1200,
                    filename = file.path(sub.dir.figures, "pub_venn_features_all.tif"), 
                    # main = "Feature sets: all data types", main.cex = 2,
                    cex = 2.8, cat.cex = 2.8, cat.pos = c(-25, 25, 180),
                    cat.dist = c(0.08, 0.08, 0.06),
                    print.mode = c("raw"), sigdigs = 1,
                    lty = "solid", lwd = 3,
                    imagetype = "tiff",
                    margin = 0.01,
                    euler.d = FALSE, scale = FALSE);

# ******* START --- create publication ready venn ***********
postscript(file = file.path(sub.dir.figures, "pub_venn_features_all.eps"), 
           onefile = FALSE, width = 6, height = 6, paper = "special", 
           family = "serif", horizontal = FALSE)
tmp <- venn.diagram(ft.set.list.total, width = 7000, height = 7000,
                   units = "px", resolution = 1200,
             filename = NULL,
             # main = "Feature sets: all data types", main.cex = 2,
             cex = 2.8, cat.cex = 2.8, cat.pos = c(-25, 25, 180),
             cat.dist = c(0.08, 0.08, 0.06),
             print.mode = c("raw"), sigdigs = 1,
             lty = "solid", lwd = 3,
             imagetype = "tiff",
             margin = 0.01,
             euler.d = FALSE, scale = FALSE);
grid.draw(tmp);
invisible(dev.off())
# ******* END ---- create publication ready venn ***********


cat("Overlap on feature level plotted in results folder.\n")

# overlap of genes (Entrez IDs) irrespective of datatype
sCCA_result_total_genes_EG = unique(unlist(sCCA_result_genes_EG))
NMF_result_total_genes_EG = unique(unlist(NMF_result_genes_EG))
MALA_result_total_genes_EG = unique(unlist(MALA_result_genes_EG))

gene.set.list.EG = list(sCCA_result_genes_EG, NMF_result_genes_EG, MALA_result_genes_EG)
names(gene.set.list.EG) = methods
gene.set.list.total.EG = list(sCCA_result_total_genes_EG, 
                              NMF_result_total_genes_EG, 
                              MALA_result_total_genes_EG)
names(gene.set.list.total.EG) = paste0(methods, " (", lapply(gene.set.list.total.EG, length), ")")

cat("Length of resulting gene sets:", names(gene.set.list.total.EG), "\n")
cat("Cardinality of the union of gene sets:", 
    length(Reduce(union, gene.set.list.total.EG)), "\n")
invisible(venn.diagram(gene.set.list.total.EG, 
             filename = file.path(sub.dir.figures, "venn_genes_all.png"),
             main = "Gene sets: all data types", main.cex = 2,
             cex = 1.6, cat.cex = 2, cat.pos = c(-25, 25, 180),
             cat.dist = c(0.06, 0.06, 0.04),
             print.mode = c("raw"), sigdigs = 1,
             lty = "solid", euler.d = FALSE, scale = FALSE))

tmp = venn.diagram(gene.set.list.total.EG, width = 7000, height = 7000,
                   units = "px", resolution = 1200, 
                   filename = file.path(sub.dir.figures, "pub_venn_genes_all.tif"), 
                   #main = "Gene sets: all data types", main.cex = 2,
                   cex = 2.8, cat.cex = 2.8, cat.pos = c(-25, 25, 180),
                   cat.dist = c(0.08, 0.08, 0.06),
                   print.mode = c("raw"), sigdigs = 1,
                   lty = "solid", lwd = 3,
                   imagetype = "tiff",
                   margin = 0.01,
                   euler.d = FALSE, scale = FALSE)

# ******* START --- create publication ready venn ***********
postscript(file = file.path(sub.dir.figures, "pub_venn_genes_all.eps"), 
           onefile = FALSE, paper = "special",
           width = 6, height = 6,
           family = "serif", horizontal = FALSE)
tmp = venn.diagram(gene.set.list.total.EG, width = 7000, height = 7000,
                   units = "px", resolution = 1200, 
             filename = NULL,
             #main = "Gene sets: all data types", main.cex = 2,
             cex = 2.8, cat.cex = 2.8, cat.pos = c(-25, 25, 180),
             cat.dist = c(0.08, 0.08, 0.06),
             print.mode = c("raw"), sigdigs = 1,
             lty = "solid", lwd = 3,
             imagetype = "tiff",
             margin = 0.01,
             euler.d = FALSE, scale = FALSE)
grid.draw(tmp)
invisible(dev.off())
# ******* END ---- create publication ready venn ***********

cat("Overlap on gene level plotted in results folder.\n")

# retrieve Gene symbols and Entrez IDs of biological datasets
datasets.all = list.files(file.path(.run.dir, "biologicalData", 
                                .current.biolog, "RData"))
datasets = datasets.all[grep("train|test", datasets.all, invert = TRUE)]
data.names = gsub(".RData", "", datasets)
data.paths = as.list(file.path(.run.dir, "biologicalData",
                               .current.biolog, "RData", datasets))
biolog.data = lapply(data.paths, readRDS)
names(biolog.data) = data.names

# get universe on feature and gene level
data.universe.ft = lapply(biolog.data, function(x) unique(rownames(x)))

# Assess statistical significance of overlap
p.val = assessSignOverlap(ft.set.list, gene.set.list.total.EG, data.universe.ft)

cat("Assess statistical significance of overlap between method results.\n",
    "Hypothesis: The overlap is significantly greater than for random",
    "gene sets of the same size. \np-value:", p.val, "\n")

#####------------------------------------------------------------------
# comparison of associated GO terms 
#####------------------------------------------------------------------

cat("Compare results on Gene Ontology level:\n")

cat("Creating parameter objects for HyperGTest ...\n")
cat("Running over-representation analysis ...\n")

# Determine gene universe 
data.universe.EG = mapply(mapFeaturesToEG, data.universe.ft, feature_des,
                          ft.col = c("EntrezID", "REF"))
cat("Number of features lacking a gene idnetifier in the universe:",
    paste(names(data.universe.ft), unlist(lapply(data.universe.EG, function(x)
      length(which(is.na(x)))))), "\n")

data.universe.EG = lapply(data.universe.EG, function(x) unique(x[which(!is.na(x))]))
cat("Size of universe in terms of Entrez IDs:", 
    paste(names(data.universe.EG), unlist(lapply(data.universe.EG, length))), "\n")

data.all.universe.EG = unique(unlist(data.universe.EG))
cat("Size of total gene universe in terms of Entrez IDs:", length(data.all.universe.EG), "\n")

# create parameter object and do GO analysis
param.conditional = FALSE
param.pval = 0.05
  
for(i in c("BP", "MF", "CC")) {
  for(j in 1:length(.data.to.integrate)){

    sets.list.EG = lapply(gene.set.list.EG, "[[", j)
    sets.list.EG = sets.list.EG[sapply(sets.list.EG, length) > 0] # remove empty sets 
  
    # create parameter objects for GO analysis 
    params.list = lapply(sets.list.EG, createParams,
                            universe = data.universe.EG[[j]], 
                            ontology = i, param.conditional,
                            param.pval)
    
    # do GO analysis for current data type
    over.list = lapply(params.list, hyperGTest)

    # save GO analysis result objects 
    saveRDS(over.list, file.path(sub.dir.RData,
                                 paste0(.data.to.integrate[j], "_GO", i, 
                                        ifelse(param.conditional, 
                                               yes = paste0("_cond_p", sub("\\.", "_", param.pval), ".RData"),
                                               no = ".RData"))))
    
  }
  
  # create parameters for GO analysis of merged gene lists 
  params.list.all = lapply(gene.set.list.total.EG, createParams, 
                          universe = data.all.universe.EG,
                          ontology = i, param.conditional, param.pval)
  
  over.list.all = lapply(params.list.all, hyperGTest)
  saveRDS(over.list.all, file.path(sub.dir.RData, 
                                   paste0("GO", i, 
                                          ifelse(param.conditional, 
                                                 yes = paste0("_result_cond_p", sub("\\.", "_", param.pval), ".RData"),
                                                 no = "_result.RData"))))

}

cat("Results of over-representation analysis saved to RData folder.\n")
cat("Plotting Venn-Diagrams ...\n")

# plot venn diagrams for each ontology
param.cat.size = 5
corr.method = "BH" # c("none", "BH")
param.FDR = 0.2
  
for(i in c("BP", "MF", "CC")){
  for(j in 1:length(.data.to.integrate)) {
    
    
    over.list = readRDS(file.path(sub.dir.RData, 
                                  paste0(.data.to.integrate[j], "_GO", i, 
                                         ifelse(param.conditional,
                                                yes = paste0("_cond_p", sub("\\.", "_", param.pval), ".RData"),
                                                no = ".RData")))) 
      
    # make number of associated GO terms accessable within result object     
    all.summary = lapply(over.list, summary, pvalue = 2, 
                         categorySize  = param.cat.size)
    
    # consider only terms with Count != 0 for the correction
    all.summary.count = lapply(all.summary, function(a) {
      subset(a, subset = a$Count != 0)
    })
    
    # do p-value correction and remove not-significant terms
    over.summary = lapply(all.summary.count, function(a) {
      p.tmp = p.adjust(a$Pvalue, corr.method, nrow(a))
      a[which(p.tmp <= param.FDR),] })
    
    # extract GO IDs of different domains
    GO.sets = lapply(over.summary, "[[", paste0("GO", i, "ID"))
    GO.sets = GO.sets[methods]
    empty.set = which(sapply(GO.sets, is.null))
    if(length(empty.set) > 0) {
      GO.sets[[empty.set]] = character(0)
    }
    
    # plot Venn diagrams for each Ontology and each data type
    names(GO.sets) = paste0(methods, " (", lapply(GO.sets, length), ")") 
    cat("Cardinality of the union of GO term for", i, "in dataset", 
        .data.to.integrate[j], ":", length(Reduce(union, GO.sets)), "\n")
    tmp = venn.diagram(GO.sets, filename = file.path(sub.dir.figures,
                                                     paste0("venn_", .data.to.integrate[j], "_GO", i, 
                                                            ifelse(param.conditional,
                                                                   yes = paste0("_cond_p", sub("\\.", "_", param.pval)),
                                                                   no = ""), 
                                                            "_catsize_", param.cat.size, 
                                                            "_corr_", corr.method, 
                                                            "_", sub("\\.", "_", param.FDR), ".png")),
                       main = paste("GO", i, ":", .data.to.integrate[j]),
                       main.cex = 2, cex = 1.6,
                       cat.cex = 2, 
                       cat.pos = c(-30,30,180), 
                       cat.dist = c(0.07, 0.07, 0.04), # alpha = c(0.35,0.35,0.5),
                       print.mode = c("raw"), sigdigs = 1,
                       lty = "solid", euler.d = FALSE, scale = FALSE)

  }
  
  over.list.all = readRDS(file.path(sub.dir.RData, 
                                    paste0("GO", i, "_result",
                                           ifelse(param.conditional,
                                                  yes = paste0("_cond_p", sub("\\.", "_", param.pval)),
                                                  no = ""), ".RData")))
  
  # make number of associated GO terms accessable within result object
  all.summary.all = lapply(over.list.all, summary, pvalue = 2,
                           categorySize = param.cat.size)
  
  # consider only terms with Count != 0 for the correction
  all.summary.all.count = lapply(all.summary.all, function(a) {
    subset(a, subset = a$Count != 0)
  })
  
  # do p-value correction and remove not-significant terms
  over.summary.all = lapply(all.summary.all.count, function(a) {
    a$Pvalue = p.adjust(a$Pvalue, corr.method, nrow(a))
    a[which(a$Pvalue <= param.FDR),] })
  
  # extract GO IDs of differnt domains
  GO.sets.all = lapply(over.summary.all, "[[", paste0("GO", i, "ID"))
  
  # plot Venn diagrams for each Ontology of merged results
  names(GO.sets.all) = paste0(methods, " (", lapply(GO.sets.all, length), ")") 
  cat("Cardinality of the union of GO term for", i, "in all datasets:", 
      length(Reduce(union, GO.sets.all)), "\n")
  tmp = venn.diagram(GO.sets.all, filename = file.path(sub.dir.figures,
                                                       paste0("venn_all_GO", i, 
                                                              ifelse(param.conditional,
                                                                     yes = paste0("_cond_p", sub("\\.", "_", param.pval)),
                                                                     no = ""),
                                                              "_catsize_", param.cat.size,
                                                              "_corr_", corr.method, 
                                                              "_", sub("\\.", "_", param.FDR), ".png")),
                     main = paste("GO", i, ": all data types"),
                     main.cex = 2, cex = 1.6,
                     cat.cex = 2, 
                     cat.pos = c(-30,30,180), 
                     cat.dist = c(0.07, 0.07, 0.04), # alpha = c(0.35,0.35,0.5),
                     print.mode = c("raw"), sigdigs = 1,
                     lty = "solid", euler.d = FALSE, scale = FALSE)
  
  tmp = venn.diagram(GO.sets.all, filename = file.path(sub.dir.figures,
                                                       paste0("pub_venn_all_GO", i,
                                                              ifelse(param.conditional,
                                                                     yes = paste0("_cond_p", sub("\\.", "_", param.pval)),
                                                                     no = ""),
                                                              "_catsize_", param.cat.size,
                                                              "_corr_", corr.method, 
                                                              "_", sub("\\.", "_", param.FDR), ".tif")), #NULL,
                     width = 7000, height = 7000, resolution = 1200,
                     #main = paste("GO", i, ": all data types"), main.cex = 2,
                     cex = 2.8, cat.cex = 2.8,
                     cat.pos = c(-25,25,180),
                     cat.dist = c(0.08, 0.08, 0.06), # alpha = c(0.35,0.35,0.5),
                     print.mode = c("raw"), sigdigs = 1,
                     lty = "solid", lwd = 3,
                     imagetype = "tiff",
                     margin = 0.01,
                     euler.d = FALSE, scale = FALSE)
  
  # ******* START --- create publication ready venn ***********
  names(GO.sets.all) = paste0(methods, " (", lapply(GO.sets.all, length), ")")
  postscript(file = file.path(sub.dir.figures,
                              paste0("pub_venn_all_GO", i,
                                     ifelse(param.conditional,
                                            yes = paste0("_cond_p", sub("\\.", "_", param.pval)),
                                            no = ""),
                                     "_catsize_", param.cat.size,
                                     "_corr_", corr.method, 
                                     "_", sub("\\.", "_", param.FDR), ".eps")),
             onefile = FALSE, paper = "special",
             width = 6, height = 6,
             family = "serif", horizontal = FALSE)
  tmp = venn.diagram(GO.sets.all, filename = NULL,
                     width = 7000, height = 7000, resolution = 1200,
                     #main = paste("GO", i, ": all data types"), main.cex = 2,
                     cex = 2.8, cat.cex = 2.8,
                     cat.pos = c(-25,25,180),
                     cat.dist = c(0.08, 0.08, 0.06), # alpha = c(0.35,0.35,0.5),
                     print.mode = c("raw"), sigdigs = 1,
                     lty = "solid", lwd = 3,
                     imagetype = "tiff",
                     margin = 0.01,
                     euler.d = FALSE, scale = FALSE)
  grid.draw(tmp)
  invisible(dev.off())
  # ******* END ----- create publication ready venn ***********

}

#####------------------------------------------------------------------
# comparison on pathway level: Topological Pathway Analysis  
#####------------------------------------------------------------------
set.seed(.biolog.seed)
cat("Compare results on Pathway level:\n")

db = "reactome"
pws = pathways("hsapiens", db)

cat("Get all pathways in", toupper(db), "(", length(pws), ") ...\n")
pw.set.name = "pw4SPIA"

cat("Prepare pathways for Topological Pathway Analysis with SPIA ...\n")
setwd(sub.dir.RData)
prepareSPIA(pws, pathwaySetName = pw.set.name)

# indicate wheather the diff should be used instead of logFC 
# (depends on data type; set to TRUE for Methylation data)
helper = rep(FALSE, length(.data.to.integrate))
helper[2] = TRUE

cat("\n", paste("-----", Sys.time(), "-----"), "\n\n")

subsets.lfc.EG.list = list()
for(i in 1:length(.data.to.integrate)){
  
  sets.list = list(sCCA_result[[i]], NMF_result[[i]], MALA_result[[i]])
  names(sets.list) = methods
  
  # calculate log ratios for result subset
  subsets.list = lapply(sets.list, function(x) {
    biolog.data[[i]][which(rownames(biolog.data[[i]]) %in% x),] })
  subsets.lfc = lapply(subsets.list, getLFC, groups = c("tumor", "normal"),
                       use.diff = helper[i])
   
  # add gene annotation information to lfc
  subsets.lfc.genes = lapply(subsets.lfc, function(lfc) {
    merge(as.data.frame(feature_des[[i]][,c("GeneSymbol", unlist(list("EntrezID", c("EntrezID", "REF"))[i]))], 
                             stringsAsFactors = FALSE), lfc,
          by.x = c("EntrezID", "REF")[i], by.y = 0, all = FALSE, sort = FALSE)
  })

  subsets.lfc.unique = lapply(subsets.lfc.genes, function(lfc) {
    lfc.order = lfc[order(abs(lfc$y), decreasing = TRUE),]
    lfc.unique = lfc.order[c(!duplicated(lfc.order$EntrezID) & !is.na(lfc.order$EntrezID)),]
    lfc.unique
  })

  # store feature lists resulting from methods with annotation
  invisible(mapply(function(x, name) {
    write.table(x, file = file.path(sub.dir.files, paste0(name, "_", 
                .data.to.integrate[i], "_result_anno.tab")),
                quote = FALSE, row.names = FALSE, sep = "\t") 
  }, subsets.lfc.genes, methods))
  
  mapply(function(x, name) {
    saveRDS(x, file = file.path(sub.dir.RData, paste0(name, "_",
            .data.to.integrate[i], "_result_anno.RData")))
  }, subsets.lfc.genes, methods)

  subsets.lfc.EG = lapply(subsets.lfc.unique, function(lfc) {
    lfc.vc = lfc$y; names(lfc.vc) = lfc$EntrezID; lfc.vc })
  
  subsets.lfc.EG.list[[i]] = subsets.lfc.EG

  # run Singnaling Pathways Impact Analysis (SPIA)
  cat("\nRun SPIA for", .data.to.integrate[i], "features and save results as RData.\n")
  res.list = mapply(function(de, all, pathwaySetName) {
    if(length(de) > 0){
      runSPIA(de, all, pathwaySetName)
    } else  {
      data.frame(Name = character(0), 
                 pSize = character(0), 
                 NDE = character(0),
                 pNDE = logical(0),
                 tA = logical(0),
                 pPERT = logical(0),
                 pG = logical(0),
                 pGFdr = numeric(0),
                 pGFWER = numeric(0),
                 Status = logical(0), 
                 stringsAsFactors = FALSE)
    }
  }, de = subsets.lfc.EG, all = data.universe.EG[i], pathwaySetName = pw.set.name, 
  SIMPLIFY = FALSE)
  
  # save SPIA result objects 
  saveRDS(res.list, file.path(sub.dir.RData, 
                              paste0(.data.to.integrate[i], "_SPIA", ".RData")))
}

# run SPIA analysis for merged method results
sets.list.merged = do.call(Map, c(c, subsets.lfc.EG.list))
sets.merged.unique = lapply(sets.list.merged, function(x) {
  sort(x)[!duplicated(names(sort(x)))]
})   

# run Singnaling Pathways Impact Analysis (SPIA)
cat("Run SPIA for merged feature set and save results as RData.\n")
res.list = lapply(sets.merged.unique, runSPIA, 
                  all = data.all.universe.EG, pathwaySetName = pw.set.name)

# save SPIA result objects 
saveRDS(res.list, file.path(sub.dir.RData, "merged_SPIA.RData"))

cat("\n", paste("-----", Sys.time(), "-----"), "\n\n")
cat("Finished SPIA.\n")
cat("Plot and compare results of SPIA. \n")

# p-value to filter pathway list for
poi.vc = c("pNDE", "pPERT", "pG", "pGFdr")
alpha.vc = c(0.05, 0.1, 0.2)

for(poi in poi.vc) {
  for(alpha in alpha.vc){
    for(j in 1:length(.data.to.integrate)){
      
      pw.list = readRDS(file.path(sub.dir.RData, 
                                  paste0(.data.to.integrate[j], 
                                         "_SPIA", ".RData")))
      # extract pathway names 
      pw.sets = lapply(pw.list, function(l) {
        if(nrow(l) > 0){
          l$Name[which(l[poi] < alpha)] 
          } else {
          character(0)
        } })
        
      
      # plot Venn diagrams of pathways for each data type
      names(pw.sets) = paste0(methods, " (", lapply(pw.sets, length), ")")
      venn.diagram(pw.sets, filename = file.path(sub.dir.figures,
                   paste0("venn_", .data.to.integrate[j], "_PW_", poi, 
                   "_", sub(".", "_", alpha, fixed = TRUE), ".png")),
                   main = paste("PWs with", poi, "<", alpha, ":", .data.to.integrate[j]),
                   main.cex = 2, cex = 1.6, cat.cex = 2, 
                   cat.pos = c(-30,30,180),
                   cat.dist = c(0.07, 0.07, 0.04), 
                   print.mode = c("raw"), sigdigs = 1, lty = "solid",
                   euler.d = FALSE, scale = FALSE)
      
    }
    
    pw.list.all = readRDS(file.path(sub.dir.RData, "merged_SPIA.RData"))
    
    # extract pathway names
    pw.sets.all = lapply(pw.list.all, function(l) {
      l$Name[which(l[poi] < alpha)] })
    
    # plot Venn diagrams of pathways for merged results
    names(pw.sets.all) = paste0(methods, " (", lapply(pw.sets.all, length), ")")
    tmp = venn.diagram(pw.sets.all, filename = file.path(sub.dir.figures,
                 paste0("venn_all_PW_", poi, "_", 
                 sub(".", "_", alpha, fixed = T), ".png")),
                 main = paste("PWs with", poi, "<", alpha, ": all data types"),
                 main.cex = 2, cex = 1.6, cat.cex = 2,
                 cat.pos = c(-30, 30, 180),
                 cat.dist = c(0.07, 0.07, 0.04),
                 print.mode = c("raw"), sigdigs = 1, lty = "solid", 
                 euler.d = FALSE, scale = FALSE)
    
    tmp = venn.diagram(pw.sets.all, filename = file.path(sub.dir.figures,
                                                         paste0("pub_venn_all_PW_", poi, "_", 
                                                                sub(".", "_", alpha, fixed = T), ".tif")), #NULL,
                       #main = paste("PWs with", poi, "<", alpha, ": all data types"),
                       #main.cex = 2, 
                       width = 7000, height = 7000, resolution = 1200,
                       cex = 2.8, cat.cex = 2.8,
                       cat.pos = c(-25, 25, 180),
                       cat.dist = c(0.08, 0.08, 0.06),
                       print.mode = c("raw"), sigdigs = 1, 
                       lty = "solid", lwd = 3,
                       margin = 0.01, imagetype = "tiff",
                       euler.d = FALSE, scale = FALSE)
    
    # ******* START --- create publication ready venn ***********
    names(pw.sets.all) = paste0(methods, " (", lapply(pw.sets.all, length), ")")
    postscript(file = file.path(sub.dir.figures,
                                paste0("pub_venn_all_PW_", poi, "_", 
                                       sub(".", "_", alpha, fixed = T), ".eps")),
               onefile = FALSE, paper = "special",
               width = 6, height = 6,
               family = "serif", horizontal = FALSE)
    tmp = venn.diagram(pw.sets.all, filename = NULL,
                 #main = paste("PWs with", poi, "<", alpha, ": all data types"),
                 #main.cex = 2, 
                 width = 7000, height = 7000, resolution = 1200,
                 cex = 2.8, cat.cex = 2.8,
                 cat.pos = c(-25, 25, 180),
                 cat.dist = c(0.08, 0.08, 0.06),
                 print.mode = c("raw"), sigdigs = 1, 
                 lty = "solid", lwd = 3,
                 margin = 0.01, imagetype = "tiff",
                 euler.d = FALSE, scale = FALSE)
    grid.draw(tmp)
    invisible(dev.off())
    # ******* END ---- create publication ready venn ***********
  }
}

#####------------------------------------------------------------------
# create tables with selected features, genes, GO terms and pathways
#####------------------------------------------------------------------

#####------------------------------------------------------------------
# Tables with features selected by each method
cat("Create tables with features and annotation information 
    for each method and data type...\n")

# load feature lists with annotation resulting from methods
ft.set.list.anno = lapply(methods, function(x) {
  paths = file.path(sub.dir.RData, paste0(x, "_", .data.to.integrate,
                                          "_result_anno.RData"))
  data = lapply(paths, readRDS)
  names(data) = .data.to.integrate
  data
})
names(ft.set.list.anno) = methods

# get gene names corresponding to Entrez IDs and add to table
ft.set.list.names = lapply(ft.set.list.anno, function(s) {
  lapply(s, function(r) { 
    if(nrow(r) > 0) {
      gene.names = select(org.Hs.eg.db, keys = r$EntrezID, columns = "GENENAME", keytype = "ENTREZID")
    } else {
	gene.names = data.frame(ENTREZID = character(0),
					GENENAME = character(0),
					stringsAsFactors = FALSE)
    } 
    names.df = cbind(r[,grep("^y$", colnames(r), invert = TRUE)], 
                     GeneName = gene.names$GENENAME, stringsAsFactors = FALSE)
  })
})
         
# print tables
mapply(function(df, n) { mapply( function(df, n) {
  if(nrow(df) > 0){
    createTables(df[order(df$GeneSymbol), ], n)
  }}, df = df, n = paste(n, names(df), "ft", sep = "_"))}, 
    df = ft.set.list.names, n = names(ft.set.list.names))


#####------------------------------------------------------------------
# Tables with GO terms associated to each method result

cat("Creating tables of over-represented GO terms",
    "for each method (sorted by size) ...\n")
cat("GO Hyper-G-Test conditional:", param.conditional, ".\n")
cat("Category size minimum:", param.cat.size, ".\n")
cat("p-value correction:", corr.method, ".\n")
cat("p-value cutoff:", param.pval, ".\n")
cat("FDR cutoff: ", param.FDR, ".\n")

for(i in c("BP", "MF", "CC")){
  for(j in 1:length(.data.to.integrate)) {
    
    over.list = readRDS(file.path(sub.dir.RData, 
                                  paste0(.data.to.integrate[j], "_GO", i, 
                                         ifelse(param.conditional,
                                                yes = paste0("_cond_p", sub("\\.", "_", param.pval), ".RData"),
                                                no = ".RData"))))
    
    # make number of associated GO terms accessable within result object
    all.summary = lapply(over.list, summary, pvalue = 2, 
                         categorySize = param.cat.size)
    
    # do p-value correction and remove not-significant terms
    over.summary = lapply(all.summary, function(a) {
      a$Padjust = p.adjust(a$Pvalue, corr.method, nrow(a))
      a[which(a$Padjust <= param.FDR),] })
    over.summary.order = lapply(over.summary, function(x) {
      x[order(x$Size, decreasing = TRUE),] })
    over.summary.sub = lapply(over.summary.order, function(x) {
      subset(x, select = grep("GO|Term|Size|Padjust", colnames(x))) })
    
    # handel potentially empty result tables
    over.summary.sub = over.summary.sub[methods]    
    empty.tab = which(sapply(over.summary.sub, is.null))
    if(length(empty.tab) > 0) {
      over.summary.sub[[empty.tab]] = data.frame(GOID = character(0),
								 Size = numeric(0),
								 Term = character(0),
								 Padjust = numeric(0),
								 stringsAsFactors = FALSE)
    }

    # plot tables: GO terms ordered by p-value 
    mapply(function(x, m) {
      items = nrow(x) #min(nrow(x), 60)
      if(items > 0){
        createTables(df = x[1:items,], fname = paste0(m, "_", .data.to.integrate[j], "_GO", i,
                                                      ifelse(param.conditional, 
                                                             yes = paste0("_cond_p", sub("\\.", "_", param.pval)),
                                                             no = ""),
                                                      "_corr_", corr.method, 
                                                      "_", sub("\\.", "_", param.FDR)))
      }
    }, x = over.summary.sub, m = methods)
  }
  
  over.list.all = readRDS(file.path(sub.dir.RData, 
                                    paste0("GO", i, "_result",
                                           ifelse(param.conditional,
                                                  yes = paste0("_cond_p", sub("\\.", "_", param.pval)),
                                                  no = ""), ".RData")))
  
  # make number of associated GO terms accessable within result object
  all.summary.all = lapply(over.list.all, summary, pvalue = 2, 
                           categorySize = param.cat.size)
  
  # do p-value correction and remove not-significant tems
  over.summary.all = lapply(all.summary.all, function(a) {
    a$Padjust = p.adjust(a$Pvalue, corr.method, nrow(a))
    a[which(a$Padjust <= param.FDR),] })
  over.summary.all.order = lapply(over.summary.all, function(x) {
    x[order(x$Size, decreasing = TRUE),] })
  over.summary.all.sub = lapply(over.summary.all.order, function(x) {
    subset(x, select = grep("GO|Term|Size|Padjust", colnames(x))) })
  
  # plot tables: GO terms ordered by p-value 
  mapply(function(x, m) {
    items = nrow(x) 
    if(items > 0){
      createTables(df = x[1:items,], fname = paste0(m, "_GO", i,
                                                    ifelse(param.conditional, 
                                                           yes = paste0("_cond_p", sub("\\.", "_", param.pval)),
                                                           no = ""),
                                                    "_corr_", corr.method, 
                                                    "_", sub("\\.", "_", param.FDR)))
    }
  }, x = over.summary.all.sub, m = methods)
}

# Merge tables of over-represented GO terms of all methods
for(i in c("BP", "MF", "CC")){
  for(j in 1:length(.data.to.integrate)) {
    
    # create file names to read
    fname = paste0(.data.to.integrate[j], "_GO", i, 
                   ifelse(param.conditional, 
                          yes = paste0("_cond_p", sub("\\.", "_", param.pval)),
                          no = ""),
                   "_corr_", corr.method, 
                   "_", sub("\\.", "_", param.FDR), ".txt")
     
    mergeTables(fname, header = c(paste0("GO", i, "ID"), "Size", "Term", "Padjust"),
                merge_by = c(paste0("GO", i, "ID"), "Size", "Term"))
  }
  
  fname = paste0("GO", i, ifelse(param.conditional, 
                                 yes = paste0("_cond_p", sub("\\.", "_", param.pval)),
                                 no = ""),
                 "_corr_", corr.method, "_", sub("\\.", "_", param.FDR), ".txt")
  
  mergeTables(fname, header = c(paste0("GO", i, "ID"), "Size", "Term", "Padjust"),
              merge_by = c(paste0("GO", i, "ID"), "Size", "Term"))
}

#####------------------------------------------------------------------
# Tables with perturbed pathways associated to each method result.

cat("Creating tables of over-represented (alpha =", alpha, 
    ")", db, "pathways for each method (sorted by size) ...\n")

for(poi in poi.vc) {
  for(alpha in alpha.vc){
    for(j in 1:length(.data.to.integrate)) {
      
      pw.list = readRDS(file.path(sub.dir.RData,
                                  paste0(.data.to.integrate[j], "_SPIA.RData")))

      pw.sets = lapply(pw.list, function(a) {
        if(nrow(a) > 0){
          a[which(a[poi] <= alpha),] 
        } else {
          a
        } })

      pw.sets.order = lapply(pw.sets, function(x) {
        x[order(x["pSize"], decreasing = TRUE),] })
      pw.sets.sub = lapply(pw.sets.order, function(x) {
        subset(x, select = grep(paste0("Name|pSize|Status|pGFdr|", poi), colnames(x))) })
      
      # plot tables: pathways ordered by size
      mapply(function(x, m) {
        items = nrow(x) #min(nrow(x), 50)
        if(items > 0){
          x = as.data.frame(x[1:items,c("Name", "pSize", poi, "pGFdr", "Status")])
          #x[,poi] = round(x[,poi], 3)
          createTables(df = x, 
                       fname = paste0(m, "_", .data.to.integrate[j], 
                                      "_PWs_", poi, "_", 
                                      sub(".", "_", alpha, fixed = TRUE),
                                      ifelse(nrow(x) > 50, yes = "", #"_top50",
                                             no = "")))
        } }, x = pw.sets.sub, m = methods, SIMPLIFY = FALSE)
    }
    
    pw.list.all = readRDS(file.path(sub.dir.RData, "merged_SPIA.RData"))
    
    pw.sets.all = lapply(pw.list.all, function(a) {
      a[which(a[poi] <= alpha),] })
    pw.sets.all.order = lapply(pw.sets.all, function(x) {
      x[order(x["pSize"], decreasing = TRUE),] })
    pw.sets.all.sub = lapply(pw.sets.all.order, function(x) {
      subset(x, select = grep(paste0("Name|pSize|Status|pGFdr|", poi), colnames(x))) })
    
    # plot tables: pathways ordered by p-value 
    mapply(function(x, m) {
      items = nrow(x) #min(nrow(x), 60)
      if(items > 0){
        x = as.data.frame(x[1:items,c("Name", "pSize", poi, "pGFdr", "Status")])
        #x[,poi] = round(x[,poi], 3)
        createTables(df = x, 
                     fname = paste0(m, "_PWs_", poi, "_", 
                                    sub(".", "_", alpha, fixed = TRUE),
                                    ifelse(nrow(x) > 60, yes = "", #"_top60",
                                           no = "")))
      } }, x = pw.sets.all.sub, m = methods, SIMPLIFY = FALSE)
  }
}

# Merge tables of perturbed pathways of all methods
for(poi in poi.vc) {
  for(alpha in alpha.vc){
    for(j in 1:length(.data.to.integrate)) {
      
      # create file names to read
      fname = paste0(.data.to.integrate[j], 
                     "_PWs_", poi, "_", 
                     sub(".", "_", alpha, fixed = TRUE),
                     ".txt")
      
      mergeTables(fname, header = c("Name", "pSize", "pGFdr", "pGFdr.1", "Status"),
                  merge_by = c("Name", "pSize"))
    }
    
    # create file names to read
    fname = paste0("PWs_", poi, "_", 
                   sub(".", "_", alpha, fixed = TRUE), 
                   ".txt")
    
    mergeTables(fname, header = c("Name", "pSize", "pGFdr", "pGFdr.1", "Status"),
                merge_by = c("Name", "pSize"))
  }
}

#####------------------------------------------------------------------
# investigate intersection with known cancer gene sets
#####------------------------------------------------------------------
set.seed(.biolog.seed)

#### KEGG, PAM50, Sorlie500, Hu306, Bedognetti, immune, target
# size of cancer gene sets
cat("Size of gene set KEGG:", length(getCancerGenes("kegg")), "\n")
cat("Size of gene set PAM:", length(getCancerGenes("pam")), "\n")
cat("Size of gene set Sorlie:", length(getCancerGenes("sorlie")), "\n")
cat("Size of gene set Hu:", length(getCancerGenes("hu")), "\n")
cat("Size of gene set Bedognetti:", length(getCancerGenes("bedo")), "\n")
cat("Size of gene set Immune:", length(getCancerGenes("immune")), "\n")
cat("Size of gene set TARGET:", length(getCancerGenes("target")), "\n")

# build intersection of cancer genesets with the gene universe
in.kegg.univ = intersect(getCancerGenes("kegg"), unlist(data.universe.EG))
in.pam50.univ = intersect(getCancerGenes("pam"), unlist(data.universe.EG))
in.sorlie500.univ = intersect(getCancerGenes("sorlie"), unlist(data.universe.EG))
in.hu306.univ = intersect(getCancerGenes("hu"), unlist(data.universe.EG))
in.bedo.univ = intersect(getCancerGenes("bedo"), unlist(data.universe.EG))
in.immune.univ = intersect(getCancerGenes("immune"), unlist(data.universe.EG))
in.target.univ = intersect(getCancerGenes("target"), unlist(data.universe.EG))

# print summary
invisible(mapply(function(name, in.univ) {
  cat("Intersection of universe with geneset", name, ":", length(in.univ), "\n") 
}, name = c("KEGG", "PAM50", "Sorlie500", "Hu306", "Bedognetti", "Immune", "Target"), 
in.univ = list(in.kegg.univ, in.pam50.univ, in.sorlie500.univ, in.hu306.univ, 
               in.bedo.univ, in.immune.univ, in.target.univ)))

# build intersection of method results with cancer genesets
in.kegg = lapply(gene.set.list.total.EG, intersect, in.kegg.univ)
in.pam50 = lapply(gene.set.list.total.EG, intersect, in.pam50.univ)
in.sorlie500 = lapply(gene.set.list.total.EG, intersect, in.sorlie500.univ)
in.hu306 = lapply(gene.set.list.total.EG, intersect, in.hu306.univ)
in.bedo = lapply(gene.set.list.total.EG, intersect, in.bedo.univ)
in.immune = lapply(gene.set.list.total.EG, intersect, in.immune.univ)
in.target = lapply(gene.set.list.total.EG, intersect, in.target.univ)

# print summary
invisible(mapply(function(name, intersection) {
  cat("Intersection of results with geneset:", name, "\n", 
    names(intersection), ":\n", 
    paste(lapply(intersection, length), collapse = "   "), "\n")
}, name = c("KEGG", "PAM50", "Sorlie500", "Hu306", "Bedognetti", "Immune", "Target"), 
   intersection = list(in.kegg, in.pam50, in.sorlie500, in.hu306, in.bedo, in.immune, in.target)))

# assess statistical significance of intersection with cancer genes with 
# permutation test and plot histogram
invisible(dev.new())
p.kegg = mapply(assessSign, gene.set.list.total.EG,
                cancer.set = "kegg",
                name = paste0("kegg_", methods), 
                overlap = in.kegg,
                universe = list(data.all.universe.EG))

p.pam = mapply(assessSign, gene.set.list.total.EG, 
               cancer.set = "pam",
               name = paste0("pam_", methods), 
               overlap = in.pam50,
               universe = list(data.all.universe.EG))

p.sorlie = mapply(assessSign, gene.set.list.total.EG, 
                  cancer.set = "sorlie",
                  name = paste0("sorlie_", methods), 
                  overlap = in.sorlie500,
                  universe = list(data.all.universe.EG))

p.hu = mapply(assessSign, gene.set.list.total.EG, 
              cancer.set = "hu",
              name = paste0("hu_", methods), 
              overlap = in.hu306,
              universe = list(data.all.universe.EG))

p.bedo = mapply(assessSign, gene.set.list.total.EG,
                cancer.set = "bedo",
                name = paste0("bedo_", methods),
                overlap = in.bedo,
                universe = list(data.all.universe.EG))

p.immune = mapply(assessSign, gene.set.list.total.EG,
                  cancer.set = "immune",
                  name = paste0("immune_", methods),
                  overlap = in.immune,
                  universe = list(data.all.universe.EG))

p.target = mapply(assessSign, gene.set.list.total.EG,
                  cancer.set = "target",
                  name = paste0("target_", methods),
                  overlap = in.target,
                  universe = list(data.all.universe.EG))


# print summary
invisible(mapply(function(name, p.val) {
  cat("Assess statistical significance of overlap between method results",
    "and cancer genes in:", name, "\nHypothesis: The overlap is significantly",
    "higher than for a random gene set. \np-values:", p.val, "\n") 
  }, name = c("KEGG", "PAM50", "Sorlie500", "Hu306", "Bedognetti", "Immune", "Target"), 
     p.val = list(p.kegg, p.pam, p.sorlie, p.hu, p.bedo, p.immune, p.target)))

# assess over-representation with Fisher's exact test
fp.kegg = mapply(assessOdds, gene.set.list.total.EG,
                cancer.set = "kegg",
                universe = list(data.all.universe.EG))

fp.pam = mapply(assessOdds, gene.set.list.total.EG, 
               cancer.set = "pam",
               universe = list(data.all.universe.EG))

fp.sorlie = mapply(assessOdds, gene.set.list.total.EG, 
                  cancer.set = "sorlie",
                  universe = list(data.all.universe.EG))

fp.hu = mapply(assessOdds, gene.set.list.total.EG, 
              cancer.set = "hu",
              universe = list(data.all.universe.EG))

fp.bedo = mapply(assessOdds, gene.set.list.total.EG,
                cancer.set = "bedo",
                universe = list(data.all.universe.EG))

fp.immune = mapply(assessOdds, gene.set.list.total.EG,
                  cancer.set = "immune",
                  universe = list(data.all.universe.EG))

fp.target = mapply(assessOdds, gene.set.list.total.EG,
                   cancer.set = "target",
                   universe = list(data.all.universe.EG))

# print summary
invisible(mapply(function(name, p.val) {
  cat("Assess statistical significance of Odds ratio", 
      "(over-representation) \nof method results",
      "and cancer genes in:", name, 
      "\nHypothesis: True Odds ratio is greater than 1.", p.val, "\n") 
}, name = c("KEGG", "PAM50", "Sorlie500", "Hu306", "Bedognetti", "Immune", "Target"), 
p.val = list(fp.kegg, fp.pam, fp.sorlie, fp.hu, fp.bedo, fp.immune, fp.target)))

#####------------------------------------------------------------------
# get cancer signatures Entrez IDs for multi-venn diagrams

# draw Venn diagrams for each method result with selected cancer signatures
sCCA.signatures = list(gene.set.list.total.EG$sCCA, in.kegg.univ, 
                       in.pam50.univ, in.target.univ)
names(sCCA.signatures) = c("sCCA", "KEGG", "PAM", "Target")

NMF.signatures = list(gene.set.list.total.EG$NMF, in.kegg.univ, 
                      in.pam50.univ, in.target.univ)
names(NMF.signatures) = c("NMF", "KEGG", "PAM", "Target")

MALA.signatures = list(gene.set.list.total.EG$MALA, in.kegg.univ, 
                       in.pam50.univ, in.target.univ)
names(MALA.signatures) = c("MALA", "KEGG", "PAM", "Target")

tmp = venn.diagram(sCCA.signatures,
             filename = file.path(sub.dir.figures, paste0("venn_sCCA_signatures.png")),
             cex = 1.6, cat.cex = 2, col = c("red", "black", "black", "black"),
             cat.col = c("red", rep("black", 3)), #cat.dist = c(0.07, 0.07, 0.05),
             label.col = c("black", "black", "black", "red", "red", "red", "black", 
                           "black", "red", "red", "red", "red", "black", 
                           "black", "red"),
             print.mode = c("raw"), sigdigs = 1,
             lty = "solid") #cat.pos = c(-30,30,180))

tmp = venn.diagram(NMF.signatures,
             filename = file.path(sub.dir.figures, paste0("venn_NMF_signatures.png")),
             cex = 1.6, cat.cex = 2, col = c("red", "black", "black", "black"),
             cat.col = c("red", rep("black", 3)), #cat.dist = c(0.07, 0.07, 0.05),
             label.col = c("black", "black", "black", "red", "red", "red", "black", 
                           "black", "red", "red", "red", "red", "black", 
                           "black", "red"),
             print.mode = c("raw"), sigdigs = 1,
             lty = "solid") #cat.pos = c(-30,30,180))

tmp = venn.diagram(MALA.signatures,
             filename = file.path(sub.dir.figures, paste0("venn_MALA_signatures.png")),
             cex = 1.6, cat.cex = 2, col = c("red", "black", "black", "black"),
             cat.col = c("red", rep("black", 3)), #cat.dist = c(0.07, 0.07, 0.05),
             label.col = c("black", "black", "black", "red", "red", "red", "black", 
                           "black", "red", "red", "red", "red", "black", 
                           "black", "red"),
             print.mode = c("raw"), sigdigs = 1,
             lty = "solid") #cat.pos = c(-30,30,180))

# venn diagrams of selected cancer signatures with all three methods
kegg.signatures = c(gene.set.list.total.EG, list(in.kegg.univ))
names(kegg.signatures) = c("sCCA (934)    ",  "  NMF (951)",  "MALA (387)  ",  "  KEGG (306)")

pam.signatures = c(gene.set.list.total.EG, list(in.pam50.univ))
names(pam.signatures) = c("sCCA (934)    ",  "  NMF (951)",  "MALA (387)  ",  "  PAM (49)")

target.signatures = c(gene.set.list.total.EG, list(in.target.univ))
names(target.signatures) = c("sCCA (934)    ",  "  NMF (951)",  "MALA (387)  ",  " Target (135)")

tmp = venn.diagram(kegg.signatures,
             filename = file.path(sub.dir.figures, paste0("venn_KEGG_signatures.png")),
             cex = 1.6, cat.cex = 1.8, col = c("black", "black", "red", "black"),
             cat.col = c(rep("black", 3), "red"), #cat.dist = c(0.1, 0.1, 0.1, 0.1),
             #cat.pos = c(-23, -15, 45, 30),
             label.col = c("black", "red", "red", "black", "red", "red", "red", "red",
                           "black", "red", "red", "black", "black", "black", "black"),
             print.mode = c("raw"), sigdigs = 1, margin = 0.1,
             lty = "solid")

# ******* START --- create publication ready venn ***********
tmp <- venn.diagram(kegg.signatures, width = 10000, height = 9000,
                    units = "px", resolution = 1200,
                    filename = file.path(sub.dir.figures, "venn_KEGG_signatures.tif"), 
                    # main = "Feature sets: all data types", main.cex = 2,
                    cex = 2.4, cat.cex = 2.4, 
                    col = c("black", "black", "red", "black"),
                    cat.col = c(rep("black", 3), "red"), 
                    #cat.dist = c(0.1, 0.1, 0.1, 0.1),
                    #cat.pos = c(-23, -15, 45, 30),
                    label.col = c("black", "red", "red", "black", "red", "red", "red", "red",
                                  "black", "red", "red", "black", "black", "black", "black"),
                    print.mode = c("raw"), sigdigs = 1,
                    lty = "solid", lwd = 3,
                    imagetype = "tiff",
                    margin = 0.1,
                    euler.d = FALSE, scale = FALSE);


postscript(file = file.path(sub.dir.figures, "venn_KEGG_signatures.eps"), 
           onefile = FALSE, width = 10, height = 9, paper = "special", 
           family = "serif", horizontal = FALSE)
tmp <- venn.diagram(kegg.signatures, width = 10000, height = 9000,
                    units = "px", resolution = 1200,
                    filename = NULL,
                    # main = "Feature sets: all data types", main.cex = 2,
                    cex = 2.8, cat.cex = 2.8, 
                    col = c("black", "black", "red", "black"),
                    cat.col = c(rep("black", 3), "red"), 
                    #cat.dist = c(0.1, 0.1, 0.1, 0.1),
                    #cat.pos = c(-23, -15, 45, 30),
                    label.col = c("black", "red", "red", "black", "red", "red", "red", "red",
                                  "black", "red", "red", "black", "black", "black", "black"),
                    print.mode = c("raw"), sigdigs = 1,
                    lty = "solid", lwd = 3,
                    imagetype = "tiff",
                    margin = 0.1,
                    euler.d = FALSE, scale = FALSE);
grid.draw(tmp);
invisible(dev.off())
# ******* END --- create publication ready venn ***********

tmp = venn.diagram(pam.signatures,
             filename = file.path(sub.dir.figures, paste0("venn_PAM_signatures.png")),
             cex = 1.6, cat.cex = 1.8, 
             col = c("black", "black", "red", "black"),
             cat.col = c(rep("black", 3), "red"), 
             #cat.dist = c(0.1, 0.1, 0.1, 0.1),
             #cat.pos = c(-23, -15, 45, 30),
             label.col = c("black", "red", "red", "black", "red", "red", "red", "red",
                           "black", "red", "red", "black", "black", "black", "black"),
             print.mode = c("raw"), sigdigs = 1, margin = 0.1,
             lty = "solid")

# ******* START --- create publication ready venn ***********
tmp <- venn.diagram(pam.signatures, width = 10000, height = 9000,
                    units = "px", resolution = 1200,
                    filename = file.path(sub.dir.figures, "venn_PAM_signatures.tif"), 
                    # main = "Feature sets: all data types", main.cex = 2,
                    cex = 2.4, cat.cex = 2.4, 
                    col = c("black", "black", "red", "black"),
                    cat.col = c(rep("black", 3), "red"), 
                    #cat.dist = c(0.1, 0.1, 0.1, 0.1),
                    #cat.pos = c(-23, -15, 45, 30),
                    label.col = c("black", "red", "red", "black", "red", "red", "red", "red",
                                  "black", "red", "red", "black", "black", "black", "black"),
                    print.mode = c("raw"), sigdigs = 1,
                    lty = "solid", lwd = 3,
                    imagetype = "tiff",
                    margin = 0.1,
                    euler.d = FALSE, scale = FALSE);


postscript(file = file.path(sub.dir.figures, "venn_PAM_signatures.eps"), 
           onefile = FALSE, width = 10, height = 9, paper = "special", 
           family = "serif", horizontal = FALSE)
tmp <- venn.diagram(pam.signatures, width = 10000, height = 9000,
                    units = "px", resolution = 1200,
                    filename = NULL,
                    # main = "Feature sets: all data types", main.cex = 2,
                    cex = 2.8, cat.cex = 2.8, 
                    col = c("black", "black", "red", "black"),
                    cat.col = c(rep("black", 3), "red"), 
                    #cat.dist = c(0.1, 0.1, 0.1, 0.1),
                    #cat.pos = c(-23, -15, 45, 30),
                    label.col = c("black", "red", "red", "black", "red", "red", "red", "red",
                                  "black", "red", "red", "black", "black", "black", "black"),
                    print.mode = c("raw"), sigdigs = 1,
                    lty = "solid", lwd = 3,
                    imagetype = "tiff",
                    margin = 0.1,
                    euler.d = FALSE, scale = FALSE);
grid.draw(tmp);
invisible(dev.off())
# ******* END --- create publication ready venn ***********

tmp = venn.diagram(target.signatures,
             filename = file.path(sub.dir.figures, paste0("venn_Target_signatures.png")),
             cex = 1.6, cat.cex = 1.8, 
             col = c("black", "black", "red", "black"),
             cat.col = c(rep("black", 3), "red"), 
             #cat.dist = c(0.1, 0.1, 0.1, 0.1),
             #cat.pos = c(-23, -15, 45, 30),
             label.col = c("black", "red", "red", "black", "red", "red", "red", "red",
                           "black", "red", "red", "black", "black", "black", "black"),
             print.mode = c("raw"), sigdigs = 1, margin = 0.1,
             lty = "solid")

# ******* START --- create publication ready venn ***********
tmp <- venn.diagram(target.signatures, width = 10000, height = 9000,
                    units = "px", resolution = 1200,
                    filename = file.path(sub.dir.figures, "venn_Target_signatures.tif"), 
                    # main = "Feature sets: all data types", main.cex = 2,
                    cex = 2.4, cat.cex = 2.4, 
                    col = c("black", "black", "red", "black"),
                    cat.col = c(rep("black", 3), "red"), 
                    #cat.dist = c(0.1, 0.1, 0.1, 0.1),
                    #cat.pos = c(-23, -15, 45, 30),
                    label.col = c("black", "red", "red", "black", "red", "red", "red", "red",
                                  "black", "red", "red", "black", "black", "black", "black"),
                    print.mode = c("raw"), sigdigs = 1,
                    lty = "solid", lwd = 3,
                    imagetype = "tiff",
                    margin = 0.1,
                    euler.d = FALSE, scale = FALSE);


postscript(file = file.path(sub.dir.figures, "venn_Target_signatures.eps"), 
           onefile = FALSE, width = 10, height = 9, paper = "special", 
           family = "serif", horizontal = FALSE)
tmp <- venn.diagram(target.signatures, width = 10000, height = 9000,
                    units = "px", resolution = 1200,
                    filename = NULL,
                    # main = "Feature sets: all data types", main.cex = 2,
                    cex = 2.8, cat.cex = 2.8, 
                    col = c("black", "black", "red", "black"),
                    cat.col = c(rep("black", 3), "red"), 
                    #cat.dist = c(0.1, 0.1, 0.1, 0.1),
                    #cat.pos = c(-23, -15, 45, 30),
                    label.col = c("black", "red", "red", "black", "red", "red", "red", "red",
                                  "black", "red", "red", "black", "black", "black", "black"),
                    print.mode = c("raw"), sigdigs = 1,
                    lty = "solid", lwd = 3,
                    imagetype = "tiff",
                    margin = 0.1,
                    euler.d = FALSE, scale = FALSE);
grid.draw(tmp);
invisible(dev.off())
# ******* END --- create publication ready venn ***********

#####------------------------------------------------------------------
# Compare performace measures
#####------------------------------------------------------------------

# load TPR, FPR and ACC for methods 
perf_meas = c("TPR", "FPR", "ACC")

ACC_arr = loadPerfMeasure("ACC") * 100
TPR_arr = loadPerfMeasure("TPR") * 100 
SPC_arr = (1-loadPerfMeasure("FPR")) * 100

# create table for ggplot input
ACC_tab = melt(data.frame(Method = factor(rep(methods, each = 2), levels = methods),
                       set = factor(c("tr", "te"), levels = c("tr", "te")),
                       measure = "ACC", ACC_arr), 
               id.vars = c("Method", "set", "measure"))
TPR_tab = melt(data.frame(Method = factor(rep(methods, each = 2), levels = methods),
                       set = factor(c("tr", "te"), levels = c("tr", "te")),
                       measure = "TPR", TPR_arr),
               id.vars = c("Method", "set", "measure"))
SPC_tab = melt(data.frame(Method = factor(rep(methods, each = 2), levels = methods),
                          set = factor(c("tr", "te"), levels = c("tr", "te")),
                          measure = "SP", SPC_arr),
               id.vars = c("Method", "set", "measure"))
                    
  
perf_ggplot = do.call(rbind, list(ACC_tab, TPR_tab, SPC_tab))
perf_ggplot$set = factor(perf_ggplot$set, labels = c("Training Set", "Test Set"))

# create plot of performance measures for methods
plotPerfMeasureBP(perf_ggplot)

cat("Comparison of methods on biological data finished.\n")

sessionInfo()

sink()
close(log.con)




