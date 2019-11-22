# Updated Multi-Omics p-Values via MiniMax Statistic
# Gabriel Odom
# 2019-11-19

library(tidyverse)



######  Find the Simulation Results  ##########################################
# top_dir <- "C:/Users/godom/" # FIU
# top_dir <- "C:/Users/gjo15/" # UM
top_dir <- "/Users/gabrielodom/" # Mac
proj_dir <- "Dropbox (BBSR)/pathwayPCA_multiOmics/"
res_dir <- "pathwayPCA_multiOmics_synthetic/sim_20190708/"

resultFiles_char <- list.files(
  path = paste0(top_dir, proj_dir, res_dir)
)

# data_ls <- readRDS(
#   paste0(top_dir, proj_dir, res_dir, resultFiles_char[1600])
# )



######  Add the MiniMax  ######################################################
# data_df <- 
#   data_ls$multiOmicsSignif %>% 
#   rowwise() %>% 
#   mutate(miniMax = max(ome1, ome2, na.rm = TRUE)) 
# 
# nullP_num <-
#   data_df %>% 
#   filter(!DEpath) %>% 
#   pull(miniMax)
# 
# hist(nullP_num, xlim = c(0, 1), main = "MiniMax p-Value of Noise")
# 
# mean(nullP_num); var(nullP_num)
# beta <- 2 * (1 / mean(nullP_num) - 1)
# 
# hist(rbeta(1000, shape1 = 2, shape2 = beta))
# 
# data2_df <- 
#   data_df %>% 
#   mutate(pMiniMax = qbeta(miniMax, shape1 = 2, shape2 = beta)) %>% 
#   select(-pathway, -interactMod) %>% 
#   arrange(pMiniMax)
# data2_df

# Well, this just sucks. We might as well be random guessing, but this may be
#   worse. We can't detect anything unless *both* data sets are close to
#   significant.
# I think we need a different simulation study...

# UPDATE: come on Gabriel. Stop being so pessimistic. Just check the results
#   for a strong signal case. The miniMax is doing just fine.
# Perhaps the reason this isn't working is due to the independence of the two
#   data platforms?



######  Make a Function  ######################################################
tabulateSimResults <- function(designFile_char){
  
  data_ls <- readRDS(
    paste0(top_dir, proj_dir, res_dir, designFile_char)
  ) 
  
  
  ###  Table of p-Values  ###
  designSpecs_df <- data.frame(
    t(data_ls$parameters[c("pctDE", "pctEffectSize")])
  )
  
  data_df <- 
    data_ls$multiOmicsSignif %>% 
    rowwise() %>% 
    mutate(miniMax = max(ome1, ome2, na.rm = TRUE)) %>% 
    bind_cols(designSpecs_df[rep(1, nrow(.)), ]) %>% 
    select(pctDE, pctEffectSize, DEpath, everything()) %>% 
    select(-pathway, -interactMod)
  
  
  ###  Beta Distribution Parameters for miniMax  ###
  trueAlpha <- 2
  trueBeta <- length(unique(data_ls$classification$DataType)) - 1
  
  xBar <- 
    data_df %>% 
    filter(!DEpath) %>% 
    pull(miniMax) %>% 
    mean()
  estBeta <- 2 * (1 / xBar - 1)
  
  
  ###  Add MiniMax p-Values and FDRs  ###
  # We want to estimate the p-values with the asymptotic estimate of the 
  #   distribution as well as the small-sample estimate
  data2_df <- 
    data_df %>% 
    mutate(
      pMiniMaxAsy = qbeta(
        miniMax, shape1 = trueAlpha, shape2 = trueBeta
      )
    ) %>% 
    mutate(asyFDR = p.adjust(pMiniMaxAsy)) %>% 
    mutate(
      pMiniMaxEst = qbeta(
        miniMax, shape1 = trueAlpha, shape2 = estBeta
      )
    ) %>% 
    mutate(estFDR = p.adjust(pMiniMaxEst)) %>% 
    arrange(miniMax)
  
  
  ###  Add Decision Indicators  ###
  data2_df %>% 
    mutate(
      asyMiniMaxCorrect = case_when(
        pMiniMaxAsy < 0.05 & !DEpath ~ "TypeI",
        pMiniMaxAsy > 0.05 & DEpath ~ "TypeII",
        TRUE ~ "True"
      )
    ) %>% 
    mutate(
      estMiniMaxCorrect = case_when(
        pMiniMaxEst < 0.05 & !DEpath ~ "TypeI",
        pMiniMaxEst > 0.05 & DEpath ~ "TypeII",
        TRUE ~ "True"
      )
    )
  
}

# test
tabulateSimResults(designFile_char = resultFiles_char[800]) %>% View
tabulateSimResults(designFile_char = resultFiles_char[100]) %>% View



######  Apply the Function  ###################################################
system.time(
  pathwayPCA_miniMax_results_df <- map_dfr(
    .x = resultFiles_char,
    .f = tabulateSimResults
  )
)
# 20.678 seconds on mac

write_csv(
  pathwayPCA_miniMax_results_df,
  path = paste0(
    top_dir, proj_dir, "pathwayPCA_multiOmics_synthetic/",
    "data/pathwayPCA_results/", "miniMax_results_100runs_20191122.csv"
  )
)
# for some strange reason, write_csv either save "True" as TRUE or Excel
#   interprets "True" as TRUE. I don't know which.
# write_csv(tibble(test = "TRUE", test2 = "True", test3 = TRUE), "test.csv")
