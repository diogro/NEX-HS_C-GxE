library(GxEMM)
ldak_loc  <- "/Genomics/grid/users/damelo/.local/bin/ldak5.linux "

library(plyr)
library(doMC)
registerDoMC(64)

expr_df = t(read.table("/Genomics//ayroleslab2/lamaya//bigProject//Datasets//GXEpaper//GeneCounts/VOOMCounts_CPM1_head_hsctrl_covfree_4svs_CORRECT_Jan8.21.txt",
                       row.names = 1))
rownames(expr_df) = gsub("^X", "", rownames(expr_df))

covariates = read.table("/Genomics//ayroleslab2/lamaya//bigProject//Datasets//GXEpaper/Covariates_forGEMMA_Jan82021.txt")
ID_C = row.names(covariates[covariates$treatment == 1,])
ID_HS = row.names(covariates[covariates$treatment == 2,])

Z_HS = as.numeric(rownames(expr_df) %in% ID_HS)
Z_C = as.numeric(rownames(expr_df) %in% ID_C)

Z_all = cbind(Z_HS, Z_C)
X     <- Z_all[,-1] # fixed effect covariates. Must include Z! Column is dropped here so that cbind(1,X) is full rank

GRM = as.matrix(read.table("/Genomics//ayroleslab2/lamaya//bigProject//Datasets//GXEpaper/GRM/head.hsc.final.norelatedness1.miss50.sXX.txt"))



out_free_HSC = llply(colnames(expr_df), function(i) GxEMM(expr_df[,i], X = X, K = GRM, Z = Z_all, 
                                             gtype='free', etype='free', ldak_loc=ldak_loc),
                     .parallel = TRUE)

saveRDS(out_free_HSC, file = "/Genomics/ayroleslab2/diogro/projects/NEX-HS_C-GxE/data/output/GxEmm_all_genes_outlist.rds")
