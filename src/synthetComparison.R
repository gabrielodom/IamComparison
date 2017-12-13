#######################################################################
##### file: synthetComparison.R                                   #####
##### input: flat lists of features resulting from the            ##### 
#####        application of methods to synthetic datasets         #####
##### output: boxplots of MCC for different simulation settings   #####
##### packages: ggplot2, scales                                   #####
##### author: B. Pucher                                           #####  
##### date created: 26/04/2017                                    #####
##### last change:  10/11/2017                                    #####
#######################################################################
rm(list = ls())

# CRAN
#install.packages("ggplot2")
#install.packages("scales")

library(ggplot2)
library(scales)
library(Cairo)

#######################################################################
#####                                                             #####
##### FUNCTIONS                                                   #####   
#####                                                             #####
#######################################################################

#####------------------------------------------------------------------
# Create box plots of performance from stats runs of each parameter set.  
#####------------------------------------------------------------------
boxPlotAUC = function(sCCA.run, NMF.run, MALA.run, stats.runs, 
                      measure = "mcc", anno = "PW_") { # anno = "PW_"
  
  # parameter help:
  # measure = c("auc", "acc", "f1", "kappa", "mcc")
  
  sCCA.auc = list()
  NMF.auc = list()
  MALA.auc = list()
  
  for(run in 1:stats.runs){

    # load arrays of performance measures from method result folders
    run.sCCA.auc = readRDS(list.files(file.path(.run.dir, "synthetsCCA", 
                                                sCCA.run, "RData", paste0("run", run)),
                                      pattern = paste0("sCCA_", anno, measure, "_arr"), 
                                      full.names = TRUE))
    sCCA.auc = c(sCCA.auc, list(run.sCCA.auc))
    
    run.NMF.auc = readRDS(list.files(file.path(.run.dir, "synthetNMF", 
                                               NMF.run, "RData", paste0("run", run)),
                                     pattern = paste0("NMF_", anno, measure, "_arr"), 
                                     full.names = TRUE))
    NMF.auc = c(NMF.auc, list(run.NMF.auc))
    
    run.MALA.auc = readRDS(list.files(file.path(.run.dir, "synthetMALA", 
                                                MALA.run, "RData", paste0("run", run)),
                                      pattern = paste0("MALA_", anno, measure, "_arr"), 
                                      full.names = TRUE))
    MALA.auc = c(MALA.auc, list(run.MALA.auc))
    
  }
  
  # extract dataset parameters 
  props = factor(dimnames(sCCA.auc[[1]])[[1]], levels = dimnames(sCCA.auc[[1]])[[1]])
  added = factor(dimnames(sCCA.auc[[1]])[[2]], levels = dimnames(sCCA.auc[[1]])[[2]])
  method = factor(c(rep("sCCA", dim(sCCA.auc[[1]])[1]*dim(sCCA.auc[[1]])[2]*stats.runs), 
             rep("NMF", dim(NMF.auc[[1]])[1]*dim(NMF.auc[[1]])[2]*stats.runs), 
             rep("MALA", dim(MALA.auc[[1]])[1]*dim(MALA.auc[[1]])[2]*stats.runs)), 
             levels = c("sCCA", "NMF", "MALA"))
  effect = factor(rep(c("weak", "moderate", "strong"), each = 3),
                  levels = c("weak", "moderate", "strong"),
                  labels = c("weak (\u03B4 = 0.2)", 
                             "moderate (\u03B4 = 0.4)", 
                             "strong (\u03B4 = 0.8)"))

  # concatenate rows of auc dataframe to vector
  auc = c(do.call(c, sapply(sCCA.auc, as.numeric, simplify = FALSE)), 
          do.call(c, sapply(NMF.auc, as.numeric, simplify = FALSE)),
          do.call(c, sapply(MALA.auc, as.numeric, simplify = FALSE)))
  
  # prepare table required to create ggplot boxplot
  auc.table = data.frame(group = props, Method = method,  value = auc,
                         effect = effect)
  
  props.de = .percentage.pw * .p.groups.d + .p.random*100
  
  measure_lim = switch(measure,
                       auc = c(0, 1),
                       acc = c(0, 1),
                       f1 = c(0, 1),
                       kappa = c(0, 1),
                       mcc = c(-1, 1))
  
  bp = ggplot(auc.table, aes(x = group, y = value, fill = Method, col = Method)) +
    geom_boxplot(position = position_dodge(0.85)) +
    scale_fill_manual(values = brewer_pal(palette = "Pastel2")(3)) +
    scale_color_manual(values = brewer_pal(palette = "Set2")(3)) +
    scale_x_discrete(labels = paste0(rep(props.de, length(added)), "%")) +
    scale_y_continuous(breaks = seq(measure_lim[1], measure_lim[2], 0.5),
                       limits = c(measure_lim[1], measure_lim[2]+0.1), oob = rescale_none) +
    facet_grid(~effect, scales = "free") +
    labs(list(x = bquote('Proportion of DE features ('*p[DE]*') in datasets'), 
              y = toupper(measure))) +
    ggtitle(label = paste("Effect")) +
    theme(plot.title = element_text(size = 35), 
          panel.background = element_rect(fill = "white", colour = "black"),
          panel.grid.major.y = element_line(colour = "grey50", linetype = 2),
          panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(),
          text = element_text(size = 25), axis.text = element_text(size = 25),
          legend.position = c(0.5, 1.17), legend.background = element_rect(colour = "white"),
          legend.direction = "horizontal", legend.key = element_rect(fill = "white"),
          legend.text = element_text(size = 20))


  tiff(file = file.path(sub.dir.figures, 
                        paste0("Boxplots_", anno, toupper(measure), "_", .serial, ".tif")),
       width = 12500, height = 6000, compression = "lzw", res = 1200)
  print(bp)
  dev.off()
  
  # postscript(onefile = FALSE, width = 12.5, height = 5.5, paper = "special",
  #            family = "serif", horizontal = FALSE)
  cairo_ps(file = file.path(sub.dir.figures, 
                           paste0("Boxplots_", anno, toupper(measure), "_", .serial, ".eps")),
          width = 12.5, height = 5.5, family = "serif")
  print(bp)
  dev.off()
  
  cat("Boxplots of", toupper(measure), "saved to file.\n")
}


#######################################################################
#####                                                             #####
##### MAIN SECTION                                                #####
#####                                                             #####
#######################################################################
current.run = .current.comp
synthet.run = .current.synthet

sub.dir.name = file.path("synthetComparison", current.run)
source(file.path(.src.dir, "setSubDirPath.R"))
source(file.path(.src.dir, "drawROC.R"))

log.con = file(file.path(sub.dir.files, "synthetComparisonLog.txt"))
sink(file = log.con, type = "output")
flush(log.con)

#####------------------------------------------------------------------
# Create BoxPlots of performance in stats runs and save to file
#####------------------------------------------------------------------
sCCA.run = .sCCA.run; cat("sCCA run: ", sCCA.run)
NMF.run = .NMF.run; cat("NMF run: ", NMF.run)
MALA.run = .MALA.run; cat("MALA run: ", MALA.run)

stats.runs = .stats.runs

boxPlotAUC(sCCA.run, NMF.run, MALA.run, stats.runs, measure = "mcc")

#####------------------------------------------------------------------
# Test pathways for over-representation in method results and save
# array of performance measure for each stats run. 
#####------------------------------------------------------------------
data.dir = file.path("syntheticData", synthet.run)

for(run in 1:.stats.runs){
  cat("\nCurrent stats run: ", run, "\n\n")
  
  GE_PW_idx = readRDS(file.path(.run.dir, data.dir, "RData", 
                                paste0("run", run), "GE_PWs.RData"))
  MET_PW_idx = readRDS(file.path(.run.dir, data.dir, "RData", 
                                 paste0("run", run), "MET_PWs.RData"))
  PW_idx = list(GE_PW = GE_PW_idx, MET_PW = MET_PW_idx)
  
  # sCCA_result
  load(file.path(.run.dir, "synthetsCCA", .sCCA.run, "RData",
                 paste0("run", run), "sCCA_result.RData"))
  # NMF_result  
  load(file.path(.run.dir, "synthetNMF", .NMF.run, "RData",
                 paste0("run", run), "NMF_result.RData"))  
  # MALA_result
  load(file.path(.run.dir, "synthetMALA", .MALA.run, "RData",
                 paste0("run", run), "MALA_result.RData"))
 
  
  result_list = list(sCCA_result, NMF_result, MALA_result)
  method = c("sCCA", "NMF", "MALA")
  dest = c(file.path(.run.dir, "synthetsCCA", .sCCA.run, "RData", paste0("run", run)),
           file.path(.run.dir, "synthetNMF", .NMF.run, "RData", paste0("run", run)),
           file.path(.run.dir, "synthetMALA", .MALA.run, "RData", paste0("run", run)))
  
  invisible(mapply(calculatePwAUC, result_list, list(PW_idx), method, synthet.run, 
                   run, dest, measure = "mcc"))
  
}

#####------------------------------------------------------------------
# Create BoxPlots of performance from over-representation analysis. 
#####------------------------------------------------------------------
sCCA.run = .sCCA.run; cat("sCCA run: ", sCCA.run, "\n")
NMF.run = .NMF.run; cat("NMF run: ", NMF.run, "\n")
MALA.run = .MALA.run; cat("MALA run: ", MALA.run, "\n")

stats.runs = .stats.runs

boxPlotAUC(sCCA.run, NMF.run, MALA.run, stats.runs, measure = "mcc", 
           anno = "PW_")


sink()
close(log.con)
