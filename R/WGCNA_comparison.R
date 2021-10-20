library(plyr)
library(ggplot2)
library(corpcor)
library(evolqg)
library(tidyverse)
library(cowplot)
library(patchwork)
library(psych)

if(!require(WGCNA)){BiocManager::install("WGCNA"); library(WGCNA)}
if(!require(tictoc)){install.packages("tictoc"); library(tictoc)}



#data_folder = "/Genomics/ayroleslab2/lamaya/bigProject/eQTLcatalog/modularity/matrices"
data_folder = "data"
dir(data_folder)

expr_list = list(head = read.table(file.path(data_folder,
                                             "head_table_wcgna_fdrLevel-1e-04.csv"),
                                   row.names = 1, header = TRUE, sep = ","),
                 body = read.table(file.path(data_folder,
                                             "body_table_wcgna_fdrLevel-1e-05.csv"),
                                   row.names = 1, header = TRUE, sep = ","))

dissTOM_list_bicor <-
  list(body = 1 - TOMsimilarityFromExpr(expr_list[["body"]], networkType = "signed", power = 8, corType = "bicor"),
       head = 1 - TOMsimilarityFromExpr(expr_list[["head"]], networkType = "signed", power = 6, corType = "bicor"))

dissTOM = dissTOM_list_bicor[["body"]]
hierTOM = hclust(as.dist(dissTOM),method="average");
colorDynamicTOM = labels2colors (cutreeDynamic(hierTOM,method="tree"))
body_modules = cutreeDynamic(hierTOM,method="tree")
plotDendroAndColors(hierTOM,
                    colors=data.frame(colorDynamicTOM),
                    dendroLabels = FALSE, marAll = c(1, 8, 3, 1),
                    main = "Gene dendrogram and module colors, TOM dissimilarity")

dissTOM = dissTOM_list_bicor[["head"]]
hierTOM = hclust(as.dist(dissTOM),method="average");
colorDynamicTOM = labels2colors (cutreeDynamic(hierTOM,method="tree"))
head_modules = cutreeDynamic(hierTOM,method="tree")
plotDendroAndColors(hierTOM,
                    colors=data.frame(colorDynamicTOM),
                    dendroLabels = FALSE, marAll = c(1, 8, 3, 1),
                    main = "Gene dendrogram and module colors, TOM dissimilarity")

names(head_modules) = colnames(expr_list$head)
names(body_modules) = colnames(expr_list$body)

Head_hsbm = read_csv("data/output/SBM/clustering/head_weights-spearman_fdr-1e-04_mcmc_mode_hierarchical-SBM.csv")
Head_hsbm$X1 = NULL
Head_hsbm$tissue = "head"
Head_hsbm$WGCNA = head_modules[Head_hsbm$Gene]

body_hsbm = read_csv("data/output/SBM/clustering/body_weights-spearman_fdr-1e-05_mcmc_mode_hierarchical-SBM.csv")
body_hsbm$X1 = NULL
body_hsbm$tissue = "body"
body_hsbm$WGCNA = body_modules[body_hsbm$Gene]

WGCNA_HSBM = rbind(body_hsbm, Head_hsbm)

WGCNA_HSBM %>%
  ggplot(aes(WGCNA, Degree)) +
  geom_boxplot(aes(group = WGCNA)) +
  geom_jitter(alpha = 0.2, width = 0.2) +
  facet_wrap(~tissue, scales = "free")

plot = WGCNA_HSBM %>%
  mutate(B3 = as.factor(B3), B1 = as.factor(B1)) %>%
  ggplot(aes(B2, WGCNA, color = B3)) +
  geom_jitter(width = 0.2, height = 0.2, alpha = 0.3) + facet_wrap(~tissue, scales = "free") +
  scale_y_continuous(breaks = 0:25, labels = c("Not\nclustered", 1:25)) + scale_x_continuous(breaks = 0:25) +
  theme_cowplot() + background_grid() + labs(y = "WGCNA Modules", x = "SBM\nLevel-2 blocks") +
  scale_color_discrete(name = "SBM\nLevel-3") + theme(legend.position = "bottom")
plot


save_plot("data/output/SBM/plots/WGCNA_comparison.png", plot, base_height = 7, ncol = 2, base_asp = 1.2)
save_plot("~/Dropbox/labbio/articles/NEX_BodyHead_Control-SBM/figures//WGCNA_comparison.png", plot, base_height = 5, ncol = 2, base_asp = 1.2)


for(t in unique(WGCNA_HSBM$tissue)){
  fdr = c("body"= 5, "head"=4)
  out_dir = paste0("data/output/SBM/clustering/", t, "_weights-spearman_fdr-1e-0", fdr[t],"_mcmc_mode_hierarchical-SBM_gene-blocks/WGCNA")
  if(dir.exists(out_dir)) {
    for(i in dir(out_dir, full.names = TRUE))
      file.remove(i)
  } else
    dir.create(out_dir, recursive = T)
  WGCNA_HSBM_t = WGCNA_HSBM %>% filter(tissue == t)
  for(block in unique(WGCNA_HSBM_t$WGCNA)){
    WGCNA_HSBM_t %>%
      filter(WGCNA == block) %>%
      select(Gene) %>%
      write_csv(file.path(out_dir, paste0("wgcna_", block, ".csv")), col_names = FALSE)
  }
}
