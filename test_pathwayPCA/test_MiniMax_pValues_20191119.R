# Updated Multi-Omics p-Values via MiniMax Statistic
# Gabriel Odom
# 2019-11-19

library(tidyverse)

######  Find the Simulation Results  ##########################################
top_dir <- "C:/Users/godom/"
proj_dir <- "Dropbox (BBSR)/pathwayPCA_multiOmics/"
res_dir <- "pathwayPCA_multiOmics_synthetic/sim_20190708/"

resultFiles_char <- list.files(
  path = paste0(top_dir, proj_dir, res_dir)
)

data_ls <- readRDS(
  paste0(top_dir, proj_dir, res_dir, resultFiles_char[1])
)


######  Add the MiniMax  ######################################################
data_df <- 
  data_ls$multiOmicsSignif %>% 
  rowwise() %>% 
  mutate(miniMax = max(ome1, ome2))

nullP_num <-
  data_df %>% 
  filter(!DEpath) %>% 
  pull(miniMax)

hist(nullP_num, xlim = c(0, 1), main = "MiniMax p-Value of Noise")

mean(nullP_num); var(nullP_num)
beta <- 2 * (1 / mean(nullP_num) - 1)

hist(rbeta(1000, shape1 = 2, shape2 = beta))

data2_df <- 
  data_df %>% 
  mutate(pMiniMax = qbeta(miniMax, shape1 = 2, shape2 = beta)) %>% 
  select(-pathway, -interactMod) %>% 
  arrange(pMiniMax)
data2_df

# Well, this just sucks. We might as well be random guessing, but this may be
#   worse. We can't detect anything unless *both* data sets are close to
#   significant.
# I think we need a different simulation study...
# Perhaps the reason this isn't working is due to the independence of the two
#   data platforms?
