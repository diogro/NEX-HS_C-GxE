# BiocManager::install("org.Dm.eg.db")

BiocManager::install("hfang-bristol/XGR", dependencies=T)
library(XGR)

library("AnnotationDbi")
library("org.Dm.eg.db")
# flyGO <- xRDataLoader(RData.customised = 'org.Dm.egGOBP',
#                       RData.location = "https://github.com/hfang-bristol/RDataCentre/blob/master/dnet/1.0.7")
# saveRDS(flyGO, "org.Dm.egGOBP.2021-09-27")
flyGO = readRDS("org.Dm.egGOBP.2021-09-27")

#BiocManager::install("clusterProfiler")
library(clusterProfiler)

block_path = "data/output/SBM/clustering/head_weights-spearman_fdr-1e-05_mcmc_mode_hierarchical-SBM_gene-blocks"
level = 1
block  = 11
getGeneList = function(block_number, level, folder_path){
  file_path = file.path(folder_path, paste0("Level_", level))
  file = dir(file_path, pattern = paste0("^", block_number), full.names = T)
  gene_list = read.csv(file, header = FALSE)[,1]
  background = read.csv(file.path(folder_path, "background.csv"), header = FALSE)[,1]
  fb = list(genes = gene_list, background = background)
  en = list(genes = bitr(fb$genes, fromType="FLYBASE",
                         toType="ENTREZID", OrgDb="org.Dm.eg.db")$ENTREZID,
            background = bitr(fb$background, fromType="FLYBASE",
                              toType="ENTREZID", OrgDb="org.Dm.eg.db")$ENTREZID)
  list(fb = fb, en = en)
}


x = getGeneList(6, 2, block_path)
out = xEnricherGenes(data = x$en$genes,
               ontology.customised = flyGO,
               ontology = 'NA',
               background = x$en$background,
               verbose = TRUE,
               check.symbol.identity = TRUE,
               ontology.algorithm = "none")
xEnrichViewer(out)

enGo = enrichGO(gene          = x$en$genes,
                universe      = x$en$background,
                OrgDb         = org.Dm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(enGo)
