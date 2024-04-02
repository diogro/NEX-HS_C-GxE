pak::pkg_install(c("plyr",  "ggplot2", "corpcor", "evolqg", "tidyverse", "cowplot", "patchwork", "psych", "WGCNA", "tictoc"))
isntallibrary(plyr)
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

#out_path = "B:/Dropbox/labbio/articles/SBM_manuscript/"
out_path = "~/Dropbox/labbio/articles/SBM_manuscript/"
#plots_path = "~/Dropbox/labbio/articles/NEX_BodyHead_Control-SBM/figures/"
out_path = "data/output/SBM/sparse_wgcna"
# create output directory
if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)

expr_list = list(head = read.table(file.path(data_folder,
                                             "head_table_WGCNA_fdrLevel-1e-02.csv"),
                                   row.names = 1, header = TRUE, sep = ","),
                 body = read.table(file.path(data_folder,
                                             "body_table_WGCNA_fdrLevel-1e-03.csv"),
                                   row.names = 1, header = TRUE, sep = ","))

sparse_networks = list(head = read.table("SBM/head_matrix_trimed.csv",
                                         header = FALSE, sep = ","),
                       body = read.table("SBM/body_matrix_trimed.csv",
                                         header = FALSE, sep = ","))                                   

dissTOM_list_bicor <-
  list(body = 1 - TOMsimilarityFromExpr(expr_list[["body"]], networkType = "signed", power = 8, corType = "bicor"),
       head = 1 - TOMsimilarityFromExpr(expr_list[["head"]], networkType = "signed", power = 6, corType = "bicor"))

sparse_dissTOM_list_bicor <- list(body = dissTOM_list_bicor[["body"]],
                               head = dissTOM_list_bicor[["head"]])
sparse_dissTOM_list_bicor[["body"]][sparse_networks[["body"]] == 0] <- 1
sparse_dissTOM_list_bicor[["head"]][sparse_networks[["head"]] == 0] <- 1 

pdf(paste0(out_path, "/WGCNA_dendrogram_body.pdf"), width = 3, height = 1.8, pointsize = 6)
dissTOM = dissTOM_list_bicor[["body"]]
hierTOM = hclust(as.dist(dissTOM),method="average");
colorDynamicTOM = labels2colors (cutreeDynamic(hierTOM,method="tree"))
body_modules = cutreeDynamic(hierTOM,method="tree")
plotDendroAndColors(hierTOM,
                    colors=data.frame(colorDynamicTOM),
                    dendroLabels = FALSE, marAll = c(1, 8, 3, 1),
                    main = "")
dev.off()

pdf(paste0(out_path, "/WGCNA_dendrogram_head.pdf"), width = 3, height = 1.8, pointsize = 6)
dissTOM = dissTOM_list_bicor[["head"]]
hierTOM = hclust(as.dist(dissTOM),method="average");
colorDynamicTOM = labels2colors (cutreeDynamic(hierTOM,method="tree"))
head_modules = cutreeDynamic(hierTOM,method="tree")
plotDendroAndColors(hierTOM,
                    colors=data.frame(colorDynamicTOM),
                    dendroLabels = FALSE, marAll = c(1, 8, 3, 1),
                    main = "")
dev.off()

pdf(paste0(out_path, "/WGCNA_sparse_dendrogram_body.pdf"), width = 3, height = 1.8, pointsize = 6)
dissTOM = sparse_dissTOM_list_bicor[["body"]]
hierTOM = hclust(as.dist(dissTOM),method="average");
colorDynamicTOM = labels2colors (cutreeDynamic(hierTOM,method="tree"))
sparse_body_modules = cutreeDynamic(hierTOM,method="tree")
plotDendroAndColors(hierTOM,
                    colors=data.frame(colorDynamicTOM),
                    dendroLabels = FALSE, marAll = c(1, 8, 3, 1),
                    main = "")
dev.off()

pdf(paste0(out_path, "/WGCNA_sparse_dendrogram_head.pdf"), width = 3, height = 1.8, pointsize = 6)
dissTOM = sparse_dissTOM_list_bicor[["head"]]
hierTOM = hclust(as.dist(dissTOM),method="average");
colorDynamicTOM = labels2colors (cutreeDynamic(hierTOM,method="tree"))
sparse_head_modules = cutreeDynamic(hierTOM,method="tree")
plotDendroAndColors(hierTOM,
                    colors=data.frame(colorDynamicTOM),
                    dendroLabels = FALSE, marAll = c(1, 8, 3, 1),
                    main = "")
dev.off()


names(head_modules) = colnames(expr_list$head)
names(body_modules) = colnames(expr_list$body)

names(sparse_head_modules) = colnames(expr_list$head)
names(sparse_body_modules) = colnames(expr_list$body)

table(head_modules)[-1] |> sum()
table(body_modules)[-1] |> sum()

table(sparse_head_modules)[-1] |> sum()
table(sparse_body_modules)[-1] |> sum()

Head_hsbm = read_csv("data/output/SBM/clustering/head_weights-spearman_fdr-1e-02_mcmc_mode_hierarchical-SBM.csv")
Head_hsbm$X1 = NULL
Head_hsbm$tissue = "head"
Head_hsbm$WGCNA = head_modules[Head_hsbm$Gene]
Head_hsbm$Sparse_WGCNA = sparse_head_modules[Head_hsbm$Gene]

body_hsbm = read_csv("data/output/SBM/clustering/body_weights-spearman_fdr-1e-03_mcmc_mode_hierarchical-SBM.csv")
body_hsbm$X1 = NULL
body_hsbm$tissue = "body"
body_hsbm$WGCNA = body_modules[body_hsbm$Gene]
body_hsbm$Sparse_WGCNA = sparse_body_modules[body_hsbm$Gene]

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
save_plot("SBM_E_corr.png", E_corr_plot, base_width = 5.2, base_height = 3.5)

comparison_plot = function(df) {
  ggplot(df, aes(B3, WGCNA, color = B4)) +
  geom_jitter(width = 0.2, height = 0.2, alpha = 0.3, size = 0.2) +
  scale_y_continuous(breaks = 0:25, labels = c("Not\nclustered", 1:25)) + scale_x_continuous(breaks = 0:25) +
  theme_cowplot() + background_grid() + labs(y = "WGCNA Modules", x = "SBM\nLevel-3 blocks") +
  scale_color_discrete(name = "SBM\nLevel-4") + 
  theme(axis.text.x = element_text(size = 8),
            plot.title = element_text(size = 8), 
            axis.text.y = element_text(size = 8),
            axis.title = element_text(size = 8),
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 6),
            axis.ticks.x = element_line(size = .3),
            axis.ticks.length=unit(.07, "cm")) +
  theme(legend.position = c(0.85, 0.9), 
        legend.text = element_text(size = 6), 
        legend.background = element_rect(fill ="white"), 
        axis.title.y = element_text(vjust=-10)) + 
        guides(colour = guide_legend(override.aes = list(size=2.5)))
}

comparison_plot_WGCNA = ggplot(WGCNA_HSBM, aes(WGCNA, Sparse_WGCNA)) +
  geom_jitter(width = 0.2, height = 0.2, alpha = 0.3, size = 0.2) +
  #scale_y_continuous(breaks = 0:25, labels = c("Not\nclustered", 1:25)) + 
  #scale_x_continuous(breaks = 0:25) +
  theme_cowplot() + background_grid() + labs(y = "Sparse WGCNA Modules", x = "WGCNA Modules") +
  theme(axis.text.x = element_text(size = 8),
            plot.title = element_text(size = 8), 
            axis.text.y = element_text(size = 8),
            axis.title = element_text(size = 8),
            axis.ticks.x = element_line(size = .3),
            axis.ticks.length=unit(.07, "cm")) +
  theme(axis.title.y = element_text(vjust=-10)) + 
        guides(colour = guide_legend(override.aes = list(size=2.5)))

save_plot(file.path(out_path, "WGCNA_Sparse_comparison.png"), comparison_plot_WGCNA, base_width = 5.2, base_height = 2.5)

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

save_plot("data/output/SBM/plots/WGCNA_comparison.png", plot, base_width = 5.2, base_height = 2.5)
save_plot(file.path(out_path, "/figures//WGCNA_comparison.png"), plot, base_width = 5.2, base_height = 5.2/2)

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
