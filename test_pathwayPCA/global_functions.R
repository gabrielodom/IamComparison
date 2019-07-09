# Global Helper Functions
# Gabriel Odom
# 2019-06-10

makeGroupIdx <- function(grp_idx, grpCard_int){
  
  pwy_idx <- seq.int(
    from = (grp_idx - 1) * grpCard_int + 1,
    to = grp_idx * grpCard_int
  )
  
}

makePath <- function(grp_idx, grpCard_int, pathType = c("ge_", "met_")){
  
  pathType <- match.arg(pathType)
  
  pwy_idx <- makeGroupIdx(grp_idx, grpCard_int)
  paste0(pathType, pwy_idx)
  
}

isGroupDE <- function(grp_idx, grpCard_int, DE_logi, pctDE){
  
  pwy_idx <- makeGroupIdx(grp_idx, grpCard_int)
  out_logi <- mean(DE_logi[pwy_idx]) >= pctDE
  names(out_logi) <- paste0("path", grp_idx)
  
  out_logi
  
}

decisionTable <- function(PathpVals_df, truePaths,
                          type, pValThresh = 0.01){
  
  # TP: marking a pathway as significant when it was designed as significant
  tpCount <- PathpVals_df %>% 
    filter(rawp <= pValThresh) %>% 
    select(terms) %>% 
    transmute(TP = terms %in% truePaths) %>%
    summarise(sum(TP)) %>%
    pull
  
  # FN: marking a pathway as not significant when it was designed as significant
  fnCount <- PathpVals_df %>% 
    filter(rawp > pValThresh) %>% 
    select(terms) %>% 
    transmute(TP = terms %in% truePaths) %>%
    summarise(sum(TP)) %>%
    pull
  
  # FP: marking a pathway as significant when it was designed not significant
  fpCount <- PathpVals_df %>% 
    filter(rawp <= pValThresh) %>% 
    select(terms) %>% 
    transmute(TP = !(terms %in% truePaths)) %>%
    summarise(sum(TP)) %>%
    pull
  
  # TN: marking a pathway as not significant when it was designed not significant
  tnCount <- PathpVals_df %>% 
    filter(rawp > pValThresh) %>% 
    select(terms) %>% 
    transmute(TP = !(terms %in% truePaths)) %>%
    summarise(sum(TP)) %>%
    pull
  
  data.frame(
    DataType = type,
    pValue = pValThresh,
    TP = tpCount,
    FN = fnCount,
    FP = fpCount,
    TN = tnCount,
    stringsAsFactors = FALSE
  )
  
}

pathSignif <- function(pathway, resp, omicsOut1, omicsOut2){
  # browser()
  
  path_df <- data.frame(
    Resp = resp,
    ome1 = getPathPCLs(omicsOut1, pathway)$PCs$V1,
    ome2 = getPathPCLs(omicsOut2, pathway)$PCs$V1
  )
  
  
  # Individual p-values
  ome1_mod <- glm(Resp ~ ome1, family = binomial, data = path_df)
  ome1p <- ifelse(
    test = ome1_mod$converged,
    yes = coef(summary(ome1_mod))[-1, 4],
    no = NA
  )
  
  ome2_mod <- glm(Resp ~ ome2, family = binomial, data = path_df)
  ome2p <- ifelse(
    test = ome2_mod$converged,
    yes = coef(summary(ome2_mod))[-1, 4],
    no = NA
  )
  
  
  # global p-value
  null_mod <- glm(Resp ~ 1, family = binomial, data = path_df)
  path_mod <- glm(Resp ~ ., family = binomial, data = path_df)
  # glm(y ~ .) includes the intercept, so these models are nested.
  
  path_aov <- anova(null_mod, path_mod)
  LRpVal <- pchisq(path_aov[2, 4], df = path_aov[2, 3], lower.tail = FALSE)
  
  
  # global p-value of interaction
  pathInter_mod <- glm(Resp ~ ome1 * ome2, family = binomial, data = path_df)
  pathInter_aov <- anova(null_mod, pathInter_mod)
  LRinterpVal <- pchisq(
    pathInter_aov[2, 4],
    df = pathInter_aov[2, 3],
    lower.tail = FALSE
  )
  
  
  # Return
  pVals_df <- data.frame(
    ome1 = ome1p,
    ome2 = ome2p,
    global = LRpVal,
    interactMod = LRinterpVal
  )
  rownames(pVals_df) <- "pValue"
  pVals_df
  
}


confusion <- function(pathway, allFeatures_char, deFeatures_char){
  
  # ! in pathway
  outside_char <- setdiff(allFeatures_char, pathway)
  # ! DE
  unexpressed_char <- setdiff(allFeatures_char, deFeatures_char)
  
  ###  2x2 Counts  ###
  # Genes in pw & marked as DE
  n1 <- length(
    intersect(pathway, deFeatures_char)
  )
  # Genes !(in pw) & marked as DE
  n2 <- length(
    intersect(outside_char, deFeatures_char)
  )
  # Genes in pw & !(marked as DE)
  n3 <- length(
    intersect(pathway, unexpressed_char)
  )
  # Genes !(in pw) & !(marked as DE)
  n4 <- length(
    intersect(outside_char, unexpressed_char)
  )
  
  ###  Return  ###
  c(n1, n2, n3, n4)
  
}

pathwayFisherExact <- function(counts_int){
  # browser()
  
  # Fisher Exact Test for overexpression
  fisher.test(
    matrix(counts_int, ncol = 2, byrow = TRUE),
    alternative = "greater"
  )$p.value
  
}
