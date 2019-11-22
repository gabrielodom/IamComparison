# Inspect Multi-Omics  MiniMax Statistic p-Values
# Gabriel Odom
# 2019-11-22

library(tidyverse)

# top_dir <- "C:/Users/godom/" # FIU
# top_dir <- "C:/Users/gjo15/" # UM
top_dir <- "/Users/gabrielodom/" # Mac
proj_dir <- "Dropbox (BBSR)/pathwayPCA_multiOmics/"
res_dir <- "pathwayPCA_multiOmics_synthetic/data/pathwayPCA_results/"

results_df <- read_csv(
  paste0(
    top_dir, proj_dir, res_dir,
    "miniMax_results_100runs_20191122.csv"
  )
)

###  Test Size  ###
results_df %>% 
  filter(!DEpath) %>% 
  group_by(pctDE, pctEffectSize) %>% 
  summarise(
    asyTS_p = mean(pMiniMaxAsy < 0.05),
    asyTS_FDR = mean(asyFDR < 0.05),
    estTS_p = mean(pMiniMaxEst < 0.05),
    estTS_FDR = mean(estFDR < 0.05)
  )
# We are waaaaaaay to conservative


###  Power  ###
results_df %>% 
  filter(DEpath) %>% 
  group_by(pctDE, pctEffectSize, run) %>% 
  summarise(
    asyTS_p = mean(pMiniMaxAsy < 0.05),
    asyTS_FDR = mean(asyFDR < 0.05),
    estTS_p = mean(pMiniMaxEst < 0.05),
    estTS_FDR = mean(estFDR < 0.05)
  ) %>% 
  ungroup() %>% 
  pivot_longer(
    cols = asyTS_p:estTS_FDR,
    names_to = "pValType",
    values_to = "Power"
  ) %>% 
  ggplot() +
    aes(x = pValType, y = Power) +
    geom_boxplot() +
    facet_grid(pctDE ~ pctEffectSize)
# I think we are underpowered.

