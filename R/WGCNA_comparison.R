pak::pkg_install(c("plyr",  "ggplot2", "corpcor", "evolqg", "tidyverse", "cowplot", "patchwork", "psych", "WGCNA", "tictoc"))
library(plyr)
library(ggplot2)
library(corpcor)
library(evolqg)
library(tidyverse)
library(cowplot)
library(patchwork)
library(psych)
library(WGCNA)
library(tictoc)

#data_folder = "/Genomics/ayroleslab2/lamaya/bigProject/eQTLcatalog/modularity/matrices"
data_folder = "data/output/SBM"
dir(data_folder)

out_path = "B:/Dropbox/labbio/articles/NEX_BodyHead_Control-SBM/"
#plots_path = "~/Dropbox/labbio/articles/NEX_BodyHead_Control-SBM/figures/"

expr_list = list(head = read.table(file.path(data_folder,
                                             "head_table_WGCNA_fdrLevel-1e-02.csv"),
                                   row.names = 1, header = TRUE, sep = ","),
                 body = read.table(file.path(data_folder,
                                             "body_table_WGCNA_fdrLevel-1e-03.csv"),
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

png("test.png")
dissTOM = dissTOM_list_bicor[["head"]]
hierTOM = hclust(as.dist(dissTOM),method="average");
colorDynamicTOM = labels2colors (cutreeDynamic(hierTOM,method="tree"))
head_modules = cutreeDynamic(hierTOM,method="tree")
plotDendroAndColors(hierTOM,
                    colors=data.frame(colorDynamicTOM),
                    dendroLabels = FALSE, marAll = c(1, 8, 3, 1),
                    main = "Gene dendrogram and module colors, TOM dissimilarity")
dev.off()


names(head_modules) = colnames(expr_list$head)
names(body_modules) = colnames(expr_list$body)

table(head_modules)[-1] |> sum()
table(body_modules)[-1] |> sum()

Head_hsbm = read_csv("data/output/SBM/clustering/head_weights-spearman_fdr-1e-02_mcmc_mode_hierarchical-SBM.csv")
Head_hsbm$X1 = NULL
Head_hsbm$tissue = "head"
Head_hsbm$WGCNA = head_modules[Head_hsbm$Gene]

body_hsbm = read_csv("data/output/SBM/clustering/body_weights-spearman_fdr-1e-03_mcmc_mode_hierarchical-SBM.csv")
body_hsbm$X1 = NULL
body_hsbm$tissue = "body"
body_hsbm$WGCNA = body_modules[body_hsbm$Gene]

WGCNA_HSBM = rbind(body_hsbm, Head_hsbm)
# write_csv(WGCNA_HSBM, file.path(out_path, "SI/TableS1-gene_clustering.csv"))

degree_plot = WGCNA_HSBM %>%
  ggplot(aes(WGCNA, Degree)) +
  geom_boxplot(aes(group = WGCNA)) +
  geom_jitter(alpha = 0.2, width = 0.2) +
  facet_wrap(~tissue, scales = "free") + theme_cowplot()
save_plot("WGCNA_degree.png", degree_plot, base_height = 7, ncol = 2, base_asp = 1.2)

B1_order = WGCNA_HSBM %>%
  group_by(B1, tissue ) %>%
  summarise(avg_degree = mean(Degree)) %>%
  arrange(desc(avg_degree))

degree_plot_body = WGCNA_HSBM %>%
  filter(tissue == "body") %>%
  ggplot(aes(x=reorder(B3,Degree), Degree)) +
  geom_boxplot(aes(group = B3)) +
  geom_jitter(alpha = 0.2, width = 0.2) +
  facet_wrap(~tissue, scales = "free") + theme_cowplot() + ggtitle("body")
degree_plot_head = WGCNA_HSBM %>%
  filter(tissue != "body") %>%
  ggplot(aes(x=reorder(B3,Degree), Degree)) +
  geom_boxplot(aes(group = B3)) +
  geom_jitter(alpha = 0.2, width = 0.2) +
  facet_wrap(~tissue, scales = "free") + theme_cowplot() + ggtitle("head")
degree_plot = degree_plot_body + degree_plot_head
save_plot("SBM_degree.png", degree_plot, base_height = 7, ncol = 2, base_asp = 1.2)

E_corr_plot_body = WGCNA_HSBM %>%
  filter(tissue == "body") %>%
  ggplot(aes(x=reorder(B1,E_corr), E_corr)) +
  geom_boxplot(aes(group = B1)) +
  geom_jitter(alpha = 0.2, width = 0.2) +
  facet_wrap(~tissue, scales = "free") + theme_cowplot() + ggtitle("body")
E_corr_plot_head = WGCNA_HSBM %>%
  filter(tissue != "body") %>%
  ggplot(aes(x=reorder(B1,E_corr), E_corr)) +
  geom_boxplot(aes(group = B1)) +
  geom_jitter(alpha = 0.2, width = 0.2) +
  facet_wrap(~tissue, scales = "free") + theme_cowplot() + ggtitle("head")
E_corr_plot = E_corr_plot_body + E_corr_plot_head
save_plot("SBM_E_corr.png", E_corr_plot, base_height = 7, ncol = 2, base_asp = 1.2)

comparison_plot = function(df) {
  ggplot(df, aes(B3, WGCNA, color = B4)) +
  geom_jitter(width = 0.2, height = 0.2, alpha = 0.3) +
  scale_y_continuous(breaks = 0:25, labels = c("Not\nclustered", 1:25)) + scale_x_continuous(breaks = 0:25) +
  theme_cowplot() + background_grid() + labs(y = "WGCNA Modules", x = "SBM\nLevel-3 blocks") +
  scale_color_discrete(name = "SBM\nLevel-4") + 
  theme(legend.position = c(0.9, 0.9), 
        legend.text = element_text(size = 23), 
        legend.background =  element_rect(fill ="white"), 
        axis.title.y = element_text(vjust=-10)) + 
        guides(colour = guide_legend(override.aes = list(size=5)))
}

plot_body = WGCNA_HSBM %>%
  filter(tissue == "body") %>%
  mutate(B4 = as.factor(B4), B1 = as.factor(B1)) %>%
  comparison_plot() + ggtitle("A. Body")
plot_head = WGCNA_HSBM %>%
  filter(tissue == "head") %>%
  mutate(B4 = as.factor(B4), B1 = as.factor(B1)) %>%
  comparison_plot() + ggtitle("B. Head")
plot = plot_body + plot_head  
save_plot("test.png", plot, base_height = 7, ncol = 2, base_asp = 1.2)

save_plot("data/output/SBM/plots/WGCNA_comparison.png", plot, base_height = 7, ncol = 2, base_asp = 1.)
save_plot(file.path(out_path, "/figures//WGCNA_comparison.pdf"), plot, base_height = 5, ncol = 2, base_asp = 1.)

for(t in unique(WGCNA_HSBM$tissue)){
  fdr = c("body"= 3, "head"=2)
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
