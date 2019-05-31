# Logistic Regression using the PCs from Each Pathway
# Gabriel Odom
# 20190507


###  Test on First Pathway  ###
# source the material from "test_gene_expr_data.R" and "test_gene_methyl_data.R"
#   first

# path1_df <- data.frame(
#   Type = tGene_df$Type,
#   geneExp = getPathPCLs(gene_aespcOut, "path1")$PCs$V1,
#   geneMethyl = getPathPCLs(methyl_aespcOut, "path1")$PCs$V1
# )
# 
# null_mod  <- glm(Type ~ 1, family = binomial, data = path1_df)
# path1_mod <- glm(Type ~ ., family = binomial, data = path1_df)
# summary(path1_mod)
# anova(
#   null_mod, path1_mod
# )


###  Make it a Function  ###
pathSignif <- function(pathway, resp, omicsOut1, omicsOut2){
  
  path_df <- data.frame(
    Resp = resp,
    ome1 = getPathPCLs(omicsOut1, pathway)$PCs$V1,
    ome2 = getPathPCLs(omicsOut2, pathway)$PCs$V1
  )
  
  
  ome1_mod <- glm(Resp ~ ome1, family = binomial, data = path_df)
  ome2_mod <- glm(Resp ~ ome2, family = binomial, data = path_df)
  
  pVals_df <- data.frame(
    ome1 = coef(summary(ome1_mod))[-1, 4],
    ome2 = coef(summary(ome2_mod))[-1, 4]
  )
  rownames(pVals_df) = "pValue"
  
  
  null_mod <- glm(Resp ~ 1, family = binomial, data = path_df)
  path_mod <- glm(Resp ~ ., family = binomial, data = path_df)
  # glm(y ~ .) includes the intercept, so these models are nested.
  
  path_aov <- anova(null_mod, path_mod)
  LRpVal <- pchisq(path_aov[2, 4], df = path_aov[2, 3], lower.tail = FALSE)
  pVals_df$global <- LRpVal
  
  
  pVals_df
  
}

# Test
pathSignif(
  pathway = "path1",
  resp = getResponse(gene_Omics),
  omicsOut1 = gene_aespcOut,
  omicsOut2 = methyl_aespcOut
)

pathSignif(
  pathway = "path5",
  resp = getResponse(gene_Omics),
  omicsOut1 = gene_aespcOut,
  omicsOut2 = methyl_aespcOut
)
