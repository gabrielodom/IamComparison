#######################################################################
##### file: biologicalData.R                                      #####
##### input: results/*CurrentRun*/preprocTCGAData-2/RData/        #####
#####        subsets.RData                                        #####
##### output: profiles (*.txt and *.RData) of samples from        #####
#####         common participants                                 #####
##### packages: --                                                #####
##### author: B. Pucher                                           #####  
##### date created: 23/07/2015                                    ##### 
##### last change:  09/11/2017                                    #####
#######################################################################
rm(list=ls())

#library("impute")

#####------------------------------------------------------------------
# TCGA barcode structure labeling samples
#####------------------------------------------------------------------
barcode.parts = c("Project", "TSS", "Participant", "Sample", "Analyte",
                  "Plate", "Center")

#####------------------------------------------------------------------
# tissue type code from TCGA Wiki at National Cancer Institute
#####------------------------------------------------------------------
tumor = 1
normal = 11
metastatic = 6
control = 20

#######################################################################
#####                                                             #####
##### FUNCTIONS                                                   #####   
#####                                                             #####
#######################################################################

#####------------------------------------------------------------------
# get type of dataset
#####------------------------------------------------------------------
getDataType = function(dataset){
  return(dataset$dataType)
}

#####------------------------------------------------------------------
# split sample names of TCGA dataset into 7 blocks
# example sample name: "TCGA-A8-A07B-01A-11R-A007Z-07"
#####------------------------------------------------------------------
split_names = function(names, sep = "-"){
  blocks = matrix(unlist(strsplit(names, split = sep)), 
                  nrow = length(names), byrow = TRUE,
                  dimnames = list(NULL, barcode.parts)) 
  
  return(blocks)
}

#####------------------------------------------------------------------
# split the 4th block of sample names to access tissue type
#####------------------------------------------------------------------
sampleType = function(name.blocks){
  type = matrix(unlist(strsplit(name.blocks[,"Sample"], split = "")),
                nrow = nrow(name.blocks), byrow = TRUE)
  type.fc = factor(as.numeric(paste(type[,1], type[,2], sep = "")),
                   levels = c(tumor, metastatic, normal, control),
                   labels = c("tumor", "metastatic", "normal", "control"))
  return(type.fc)
}


#####------------------------------------------------------------------
# split the Analyte-block of sample names to get data type
#####------------------------------------------------------------------
analyteType = function(name.blocks){
  type = matrix(unlist(strsplit(name.blocks[,"Analyte"], split = "")),
                nrow = nrow(name.blocks), byrow = TRUE)[,3]
  
  return(type)
}

#####------------------------------------------------------------------
# returns the participant id of a provided sample
#####------------------------------------------------------------------
getParticipantId = function(sample){
  return(split_names(sample)[,"Participant"])
}

#####------------------------------------------------------------------
# get participants common to all subsets
#####------------------------------------------------------------------
getCommonParticipantsTable = function(subset.list){
  commons.table = numeric(0)
  for(i in 1:length(subset.list)){
    commons = lapply(subset.list, intersect, subset.list[[i]])
    commons.table = cbind(commons.table, unlist(lapply(commons, length)))
  }
  colnames(commons.table) = names(subset.list)
  return(commons.table)
}


#####------------------------------------------------------------------
# extract profiles of samples with provided barcode
#####------------------------------------------------------------------
getProfiles = function(barcodes, dataset){
  profile.idx = which(colnames(dataset$Data) %in% barcodes)
  profile.set = dataset$Data[,profile.idx]
  row.names(profile.set) = dataset$Des[,grep("EntrezID|REF", colnames(dataset$Des))]
  profiles.sorted = profile.set[,order(colnames(profile.set))]
  nb.samples = length(barcodes)
  nb.participants = length(unique(getParticipantId(barcodes)))
  if(nb.samples > nb.participants){
    redundant.participant = 
      names(which(table(getParticipantId(colnames(profiles.sorted))) > 1))
    idx.redundant = which(getParticipantId(colnames(profiles.sorted)) ==
                            redundant.participant)
    # if one participant was measured several times within the same data 
    # type, keep only one sample
    profiles.sorted = profiles.sorted[,-idx.redundant[-1]] 
  }
  
  return(profiles.sorted)
}


#####------------------------------------------------------------------
# save profiles of samples common to all subsets to files
#####------------------------------------------------------------------
saveProfiles = function(profiles, filename, save.RData = TRUE){
  write.table(profiles, 
              file.path(sub.dir.files, paste(filename, ".txt", sep = "")),
              quote = F, row.names = T, col.names = T)
  cat("Profiles saved to ", filename, ".txt \n", sep = "")
  
  if(save.RData){
    saveRDS(profiles, 
            file.path(sub.dir.RData, paste(filename, ".RData", sep = "")))
    cat("Profiles saved to ", filename, ".RData \n", sep = "")
  }

}

#######################################################################
#####                                                             #####
##### MAIN SECTION                                                #####
#####                                                             #####
#######################################################################
current.run = .current.biolog
sub.dir.name = file.path("biologicalData", current.run)
source(file.path(.src.dir, "setSubDirPath.R"))

load(file.path(.run.dir, "preprocTCGAData-2", "RData", "subsets.RData")) 

log.con = file(file.path(sub.dir.files, "biologicalDataLog.txt"))
sink(file = log.con, type = "output")
flush(log.con)

cat("Continue preprocessing of biological datasets ...\n")

# split up subset list  
gene.exp.tu = 
  subsets[[which(lapply(subsets, getDataType) == "GeneExp_tu")]] 
gene.exp.no = 
  subsets[[which(lapply(subsets, getDataType) == "GeneExp_no")]]
methyl.tu =
  subsets[[which(lapply(subsets, getDataType) == "Methylation_tu")]]
methyl.no = 
  subsets[[which(lapply(subsets, getDataType) == "Methylation_no")]]

cat("Dimension of data subsets:",
    unlist(lapply(subsets, function(s) {
 paste(getDataType(s), paste(unlist(lapply(s, dim)$Data), collapse = " "))
})), sep = "\n")

#rm(subsets)

# sample names of each subset
ge.tu.samples = dimnames(gene.exp.tu$Data)[[2]]
ge.no.samples = dimnames(gene.exp.no$Data)[[2]]
met.tu.samples = dimnames(methyl.tu$Data)[[2]]
met.no.samples = dimnames(methyl.no$Data)[[2]]

# sample names of subsets split into blocks
ge.tu.samples.blocks = split_names(ge.tu.samples)
ge.no.samples.blocks = split_names(ge.no.samples)
met.tu.samples.blocks = split_names(met.tu.samples)
met.no.samples.blocks = split_names(met.no.samples)

#####------------------------------------------------------------------
# find common pariticipants across all data types and subsets
#####------------------------------------------------------------------
# extract the ID of all participants involved
tu.part = list(ge.tu = ge.tu.samples.blocks[,"Participant"],
               met.tu = met.tu.samples.blocks[,"Participant"])
               # prot.tu = prot.tu.samples.blocks[,"Participant"])
no.part = list(ge.no = ge.no.samples.blocks[,"Participant"],
               met.no = met.no.samples.blocks[,"Participant"])

subsets.part = c(tu.part, no.part)

# create table with participants common to pairs of subsets
getCommonParticipantsTable(subsets.part)

#Nb of participants common to all tumor subsets
common.tu = Reduce(intersect, tu.part) 
#Nb of participants common to all normal subsets
common.no = Reduce(intersect, no.part) 

#Nb of participants common to all subsets 
commons.all = intersect(common.tu, common.no) 

cat("Number of Participants common to all tumor Subsets: ",
    length(common.tu), "\n",
    "Number of Participants common to all normal Subsets: ",
    length(common.no), "\n",
    "Number of Participants common to all Subsets: ",
    length(commons.all), "\n", 
    sep = "")


cat("IDs of common Participants sorted: \n")
print(getParticipantId(sort(colnames(gene.exp.tu$Data)[
  which(getParticipantId(colnames(gene.exp.tu$Data)) %in% 
          commons.all)])), quote = F)


#samples from common participants (blocks)
ge.tu.com = ge.tu.samples.blocks[which(
  ge.tu.samples.blocks[,"Participant"] %in% commons.all),]
ge.no.com = ge.no.samples.blocks[which(
  ge.no.samples.blocks[,"Participant"] %in% commons.all),]
met.tu.com = met.tu.samples.blocks[which(
  met.tu.samples.blocks[,"Participant"] %in% commons.all),]
met.no.com = met.no.samples.blocks[which(
  met.no.samples.blocks[,"Participant"] %in% commons.all),]

ge.tu.com = ge.tu.com[!duplicated(ge.tu.com[,"Participant"]),]
ge.no.com = ge.no.com[!duplicated(ge.no.com[,"Participant"]),]
met.tu.com = met.tu.com[!duplicated(met.tu.com[,"Participant"]),]
met.no.com = met.no.com[!duplicated(met.no.com[,"Participant"]),]

#samples from common participants (barcode)
ge.samples.tu.com = apply(ge.tu.com, 1, paste, collapse = "-")
ge.samples.no.com = apply(ge.no.com, 1, paste, collapse = "-")
met.samples.tu.com = apply(met.tu.com, 1, paste, collapse = "-")
met.samples.no.com = apply(met.no.com, 1, paste, collapse = "-")

#####------------------------------------------------------------------
# extract profiles of samples with provided barcode
#####------------------------------------------------------------------
ge.tu.profiles = getProfiles(ge.samples.tu.com, gene.exp.tu)
ge.no.profiles = getProfiles(ge.samples.no.com, gene.exp.no)
met.tu.profiles = getProfiles(met.samples.tu.com, methyl.tu)
met.no.profiles = getProfiles(met.samples.no.com, methyl.no)

#####------------------------------------------------------------------
# recombine tumor and normal samples of each data type
#####------------------------------------------------------------------
ge.profiles = cbind(ge.tu.profiles, ge.no.profiles)
met.profiles = cbind(met.tu.profiles, met.no.profiles)

#####------------------------------------------------------------------
# in RNAseq dataset replace zeroes with NAs
# remove cases with more than 10% NAs
#####------------------------------------------------------------------
cat("In RNASeq dataset: replace zeroes with NA\n")
ge.profiles.na = ge.profiles
ge.profiles.na[which(ge.profiles.na == 0)] = NA

NAs.th = 10
cat("Remove features (rows) with more than", NAs.th, "% NAs\n")
ge.percent.na = apply(ge.profiles.na, 1, function(r) {length(which(is.na(r)))/length(r)*100})
met.percent.na = apply(met.profiles, 1, function(r) {length(which(is.na(r)))/length(r)*100})

ge.profiles.dupl = ge.profiles.na[which(ge.percent.na < NAs.th),]
met.profiles.dupl = met.profiles[which(met.percent.na < NAs.th),]

cat("Dimensions after removing cases with >", NAs.th, "% NAs:", 
    "\n GeneExp: ", dim(ge.profiles.dupl),
    "\n Methylation: ", dim(met.profiles.dupl),
    "\n")

#####------------------------------------------------------------------
# Remove duplicated rows (features with same unique identifier 
# referring to multiple gene symbols) from the datasets.
#####------------------------------------------------------------------
cat("Remove duplicated rows (features with the same unique identifier\n",
    "referring to multiple gene symbols) from datasets.\n")

ge.profiles.filter = ge.profiles.dupl[!duplicated(rownames(ge.profiles.dupl)),]
met.profiles.filter = met.profiles.dupl[!duplicated(rownames(met.profiles.dupl)),]

cat("Dimensions after removing duplicated rows:", 
    "\n GeneExp: ", dim(ge.profiles.filter),
    "\n Methylation: ", dim(met.profiles.filter),
    "\n")

#####------------------------------------------------------------------
# save Histograms of datasets before transformation as .png and .eps
#####------------------------------------------------------------------
subsets.list <- list(ge.profiles.filter, met.profiles.filter)
names(subsets.list) <- c("RNA-seq Gene Expression Levels",
                         "DNA-Methylation Levels - Beta-values") 
# print to eps
postscript(file = file.path(sub.dir.figures, "subsets_hist.eps"), 
           onefile = FALSE, width = 8, height = 3.5, paper = "special", 
           family = "serif", horizontal = FALSE)
y = c(2*10^-5, 10)
par(mfrow = c(1,2), mar = c(2,2,6,2))
invisible(sapply(seq(1,2,1), function(x, data){
  hist(data[[x]], prob = T, breaks = 10, 
       main = names(data[x]), 
       ylim = c(0,y[x]))
  lines(density(data[[x]], na.rm = TRUE), lwd = 2, col = "darkblue")
}, data = subsets.list))
par(mfrow = c(1,1), mar = c(1,1,1,1), cex = 1)
title(main = paste("Histogram before Transformation"))
invisible(dev.off())

# print to png
png(file.path(sub.dir.figures, "subsets_hist.png"), width = 3000, 
    height = 1300, res = 300)
y = c(2*10^-5, 10)
par(mfrow = c(1,2), mar = c(2,2,6,2))
invisible(sapply(seq(1,2,1), function(x, data){
  hist(data[[x]], prob = T, breaks = 10, 
       main = names(data[x]), 
       ylim = c(0,y[x]))
  lines(density(data[[x]], na.rm = TRUE), lwd = 2, col = "darkblue")
}, data = subsets.list))
par(mfrow = c(1,1), mar = c(1,1,1,1), cex = 1)
title(main = paste("Histogram before Transformation"))
invisible(dev.off())

cat("Save Histogram before transformation and imputation ... \n")

#####------------------------------------------------------------------
# impute NA GeneExp and Methylation values 
#####------------------------------------------------------------------
# Imputation of missing values with half the lowest value in each feature
cat("Impute low-expression values in RNA-Seq data with half",
    "of the lowest value in each feature.\n")
ge.profiles.im = t(apply(ge.profiles.filter, MARGIN = 1, FUN = function(x) {
  x[which(is.na(x))] = min(x, na.rm = TRUE)/2; x}))

cat("Impute missing beta values in Methylation data with half",
    "of the lowest value in each feature.\n")
met.profiles.im = t(apply(met.profiles.filter, MARGIN = 1, FUN = function(x) {
  x[which(is.na(x))] = min(x, na.rm = TRUE)/2; x}))


#####------------------------------------------------------------------
# log2 transform gene expression data  
#####------------------------------------------------------------------
pc = 0
cat("Add pseudocount of", pc, "and log2 transform RNA-seq data.\n")
ge.profiles.tr = log2(ge.profiles.im+pc)

saveProfiles(ge.profiles.tr, "GeneExp_TCGAbarcode", save.RData = FALSE)

# Reduce sample names to "tumor_R", "normal_R", "tumor_D", ..
colnames(ge.profiles.tr) = 
  paste(sampleType(split_names(colnames(ge.profiles.tr))),
        analyteType(split_names(colnames(ge.profiles.tr))), sep = "_")

#####------------------------------------------------------------------
# transform methylation data from beta- to M-values
#####------------------------------------------------------------------
cat("Transform Methylation data from beta- to M-values ...\n")
cat("Reverse sign of M-values to mirror distribution at 0.\n")

met.profiles.tr = -1*log2(met.profiles.im/(1-met.profiles.im))

saveProfiles(met.profiles.tr, "Methyl_TCGAbarcode", save.RData = FALSE)

# Reduce sample names to "tumor_R", "normal_R", "tumor_D", ..
colnames(met.profiles.tr) = 
  paste(sampleType(split_names(colnames(met.profiles.tr))),
        analyteType(split_names(colnames(met.profiles.tr))), sep = "_")

#####------------------------------------------------------------------
# save Histograms of datasets after transformation to .png and .eps
#####------------------------------------------------------------------
subsets.tr.list <- list(ge.profiles.tr, met.profiles.tr)
names(subsets.tr.list) <- c("Log2 of Gene Expression Levels",
                            "DNA-Methylation Levels - M-values")

# print to eps
postscript(file = file.path(sub.dir.figures, "subsets_hist_tr.eps"), 
           onefile = FALSE, width = 8, height = 3.5, paper = "special", 
           family = "serif", horizontal = FALSE)
y = c(0.2,0.3)
par(mfrow = c(1,2), mar = c(2,2,6,2))
invisible(sapply(seq(1,2,1), function(x, data){
  hist(data[[x]], prob = T, breaks = 20, main = names(data[x]), ylim = c(0,y[x])) 
  lines(density(data[[x]], adjust = 3), lwd = 2, col = "darkblue")
}, data = subsets.tr.list))
par(mfrow = c(1,1), mar = c(1,1,1,1), cex = 1)
title(main = paste("Histogram after Transformation and Imputation"))
invisible(dev.off())

# print to png
png(file.path(sub.dir.figures, "subsets_hist_tr.png"), width = 3000, 
    height = 1300, res = 300)
y = c(0.2,0.3)
par(mfrow = c(1,2), mar = c(2,2,6,2))
invisible(sapply(seq(1,2,1), function(x, data){
  hist(data[[x]], prob = T, breaks = 20, main = names(data[x]), ylim = c(0,y[x])) 
  lines(density(data[[x]], adjust = 3), lwd = 2, col = "darkblue")
}, data = subsets.tr.list))
par(mfrow = c(1,1), mar = c(1,1,1,1), cex = 1)
title(main = paste("Histogram after Transformation and Imputation"))
invisible(dev.off())

cat("Save Histogram after transformation and imputation ... \n")

#####------------------------------------------------------------------
# assess differences in mean expression value per gene btw. groups 
#####------------------------------------------------------------------
cat("Assess log2 fold change of gene expression:\n")
cat("Back-transform to expression levels to calculate log fold change.\n")
ge.tu.means = rowMeans(2^ge.profiles.tr[,grep("tumor", colnames(ge.profiles.tr))])
ge.no.means = rowMeans(2^ge.profiles.tr[,grep("normal", colnames(ge.profiles.tr))])

cat("Calculate log-FC: log2(mean(tumor)/mean(normal))\n")
ge.means.lfc = log2(ge.tu.means/ge.no.means)

ge.lfc.th = 3.5

ge.sig.diff = length(which(abs(ge.means.lfc) > ge.lfc.th))

cat("Thresholds for differential gene expression:\n log2 Fold Change",
    ge.lfc.th, "\n")
cat("Number of features (%) considered as differentially expressed:\n", 
    ge.sig.diff, "(", round(ge.sig.diff/nrow(ge.profiles.tr)*100, 2), "% )\n")

# Plot and save Histogram of logFC in gene expression data 
# print to eps
postscript(file = file.path(sub.dir.figures, paste0("GeneExp_logFC", ge.lfc.th, ".eps")), 
           onefile = FALSE, width = 8, height = 4.2, paper = "special", 
           family = "serif", horizontal = FALSE)
hist(ge.means.lfc, prob = T, 
     main = "Log Fold Change of gene expression in tumor vs. normal samples",
     xlim = c(-5,5), xlab = ("LFC"),
     ylim = c(0,max(density(ge.means.lfc)$y)))
axis(side = 1, at = c(-ge.lfc.th, ge.lfc.th))
lines(density(ge.means.lfc, adjust = 3), lwd = 2, col = "darkblue")
Fn = ecdf(ge.means.lfc)
abline(v = quantile(ge.means.lfc, probs = Fn(c(-ge.lfc.th, ge.lfc.th))))
text(x = c(-ge.lfc.th, ge.lfc.th), y = c(0.6,0.6), pos = 4, 
     labels = paste0(round(Fn(c(-ge.lfc.th, ge.lfc.th))*100, 2), "%"))
invisible(dev.off())

# print to png
png(file.path(sub.dir.figures, paste0("GeneExp_logFC", ge.lfc.th, ".png")), 
              width = 3000, height = 1500, res = 300)
hist(ge.means.lfc, prob = T, 
     main = "Log Fold Change of gene expression in tumor vs. normal samples",
     xlim = c(-5,5), xlab = ("LFC"),
     ylim = c(0,max(density(ge.means.lfc)$y)))
axis(side = 1, at = c(-ge.lfc.th, ge.lfc.th))
lines(density(ge.means.lfc, adjust = 3), lwd = 2, col = "darkblue")
Fn = ecdf(ge.means.lfc)
abline(v = quantile(ge.means.lfc, probs = Fn(c(-ge.lfc.th, ge.lfc.th))))
text(x = c(-ge.lfc.th, ge.lfc.th), y = c(0.6,0.6), pos = 4, 
     labels = paste0(round(Fn(c(-ge.lfc.th, ge.lfc.th))*100, 2), "%"))
invisible(dev.off())

#####------------------------------------------------------------------
# assess differences in methylation levels per site btw. groups 
#####------------------------------------------------------------------
cat("Assess difference in means of Methylation levels:\n")
cat("Calulate means of M-values.\n")
met.tu.means = rowMeans(met.profiles.tr[,grep("tumor", colnames(met.profiles.tr))])
met.no.means = rowMeans(met.profiles.tr[,grep("normal", colnames(met.profiles.tr))])

# find differentially methylated CpG sites by assessing absolute 
# difference between M-values (Du, 2010)
cat("Calculate difference of mean M-values.\n")
met.means.diff = met.tu.means - met.no.means 

met.diff.th = 2.2

met.sig.diff = length(which(abs(met.means.diff) > met.diff.th))

cat("Thresholds for differential methylation:\nabsolute difference in mean",
    met.diff.th, "\n")
cat("Number of features (%) considered as significantly different:\n", 
    met.sig.diff, "(", round(met.sig.diff/nrow(met.profiles.tr)*100, 2), "% )\n")

# print to eps
postscript(file = file.path(sub.dir.figures, paste0("Methyl_mean_diff", met.diff.th, ".eps")), 
           onefile = FALSE, width = 8, height = 4.2, paper = "special", 
           family = "serif", horizontal = FALSE)
hist(met.means.diff, prob = T, 
     main = "Difference in mean methylation level of tumor and normal samples", 
     xlim = c(-4,4), xlab = "M-value difference",
     ylim = c(0, max(density(met.means.diff)$y)))
axis(side = 1, at = c(-met.diff.th, met.diff.th))
lines(density(met.means.diff, adjust = 3), lwd = 2, col = "darkblue")
Fn = ecdf(met.means.diff)
abline(v = quantile(met.means.diff, probs = Fn(c(-met.diff.th, met.diff.th))))
text(x = c(-met.diff.th, met.diff.th), y = c(1.5,1.5), pos = 4, 
     labels = paste0(round(Fn(c(-met.diff.th, met.diff.th))*100, 2), "%"))
invisible(dev.off())

# print to png
png(file.path(sub.dir.figures,paste0("Methyl_mean_diff", met.diff.th, ".png")), 
    width = 3000, height = 1500, res = 300)
hist(met.means.diff, prob = T, 
     main = "Difference in mean methylation level of tumor and normal samples", 
     xlim = c(-4,4), xlab = "M-value difference",
     ylim = c(0, max(density(met.means.diff)$y)))
axis(side = 1, at = c(-met.diff.th, met.diff.th))
lines(density(met.means.diff, adjust = 3), lwd = 2, col = "darkblue")
Fn = ecdf(met.means.diff)
abline(v = quantile(met.means.diff, probs = Fn(c(-met.diff.th, met.diff.th))))
text(x = c(-met.diff.th, met.diff.th), y = c(1.5,1.5), pos = 4, 
     labels = paste0(round(Fn(c(-met.diff.th, met.diff.th))*100, 2), "%"))
invisible(dev.off())

#####------------------------------------------------------------------
# save profiles of samples common to all subsets to .txt and .RData
#####------------------------------------------------------------------
saveProfiles(ge.profiles.tr, "GeneExp")
saveProfiles(met.profiles.tr, "Methyl")

#####------------------------------------------------------------------
# Prepare datasets for cross validation
#####------------------------------------------------------------------
if(.prepare.CV == TRUE){
  
  source(file.path(.src.dir, "crossValidation.R"))
  
  data.names = c("GeneExp", "Methyl")
  datasets = paste0(data.names, ".RData")
  data.paths = as.list(file.path(sub.dir.RData, datasets))
  data.raw = lapply(data.paths, readRDS)
  names(data.raw) = data.names
    
  set.seed(.biolog.seed)
  randomSplit(data.raw, subsets = .subsets.val, .pr.train, 
              samples.paired = TRUE)
  
}



sink()
close(log.con)

