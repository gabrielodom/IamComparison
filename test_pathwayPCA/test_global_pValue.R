# Logistic Regression using the PCs from Each Pathway
# Gabriel Odom
# 20190507


###  Test on First Pathway  ###
# source the material from "test_gene_expr_data.R" and "test_gene_methyl_data.R"
#   first

path1_df <- data.frame(
  Type = tGene_df$Type,
  geneExp = getPathPCLs(gene_aespcOut, "path1")$PCs$V1,
  geneMethyl = getPathPCLs(methyl_aespcOut, "path1")$PCs$V1
)

null_mod  <- glm(Type ~ 1, family = binomial, data = path1_df)
path1_mod <- glm(Type ~ ., family = binomial, data = path1_df)
summary(path1_mod)
anova(
  null_mod, path1_mod
)


###  Make it a Function  ###
pathSignif <- function(pathway, resp, omicsOut1, omicsOut2){
  
  path_df <- data.frame(
    Resp = resp,
    ome1 = getPathPCLs(omicsOut1, pathway)$PCs$V1,
    ome2 = getPathPCLs(omicsOut2, pathway)$PCs$V1
  )
  
  null_mod <- glm(Resp ~ 1, family = binomial, data = path_df)
  path_mod <- glm(Resp ~ ., family = binomial, data = path_df)
  
  path_aov <- anova(null_mod, path_mod)
  LRpVal <- pchisq(path_aov[2, 4], df = path_aov[2, 3], lower.tail = FALSE)
  
  pVals_mat <- t(coef(summary(path_mod)))[4, -1, drop = FALSE]
  pVals_df <- as.data.frame(pVals_mat, row.names = "pValue")
  pVals_df$global <- LRpVal
  
  pVals_df
  
}

# Test
pathSignif(
  pathway = "path1",
  resp = tGene_df$Type,
  omicsOut1 = gene_aespcOut,
  omicsOut2 = methyl_aespcOut
)
pathSignif(
  pathway = "path1",
  resp = getResponse(gene_Omics),
  omicsOut1 = gene_aespcOut,
  omicsOut2 = methyl_aespcOut
)
