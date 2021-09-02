if(!require(plyr)){install.packages("plyr"); library(plyr)}
if(!require(tidyverse)){install.packages("tidyverse"); library(tidyverse)}
if(!require(cowplot)){install.packages("cowplot"); library(cowplot)}
if(!require(ggthemes)){install.packages("ggthemes"); library(ggthemes)}
if(!require(xkcd)){install.packages("xkcd"); library(xkcd)}
if(!require(extrafont)){install.packages("extrafont"); library(extrafont)}
if(!require(superheat)){install.packages("superheat"); library(superheat)}
if(!require(mcclust)){install.packages("mcclust"); library(mcclust)}
if(!require(patchwork)){install.packages("patchwork"); library(patchwork)}

block_df = read_csv("data/output/SBM/clustering/head_mcmc_cutoff-spearman_val-0.5_hierarchical-SBM.csv")

block_df = block_df %>%
  arrange(B5, B4, B3, B2, B1, Degree)

e_mats = vector("list", 5)
for (i in seq_along(e_mats)){
  e_matrix = read_csv(paste0("data/output/SBM/clustering/head_mcmc_cutoff-spearman_val-0.5_hierarchical-SBM_e_matrix_level",i-1,".csv"))
  e_matrix$X1 = NULL
  e_matrix = as.matrix(e_matrix)
  rownames(e_matrix) = colnames(e_matrix) = as.character(colnames(e_matrix))
  blocks = as.character(unique(as.data.frame(block_df)[,paste0("B", i)]))
  e_matrix= e_matrix[blocks, blocks]
  #e_matrix = e_matrix[-61, -61]
  diag(e_matrix) = diag(e_matrix)/2
  e_mats[[i]] = e_matrix
}


makePlot = function(e_matrix){
  colnames(e_matrix) = rownames(e_matrix) = paste0("b", rownames(e_matrix))
  melted_cormat <- reshape2::melt(e_matrix)
  melted_cormat[melted_cormat == 0] = NA
  #plot heatmap
  ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile()+
    scale_fill_viridis_c(alpha = 1)  + theme_tufte() + labs(y = "Blocks", x = "Blocks") +
    theme(axis.text.x = element_text(angle=45, vjust=0.6), legend.position = "none") 
}
plot_list = lapply(e_mats, makePlot)

all_plots = plot_grid(plotlist = plot_list)
save_plot("data/output/SBM/plots/E_matrices.png", all_plots, base_height = 5, base_asp = 1.2, ncol = 3, nrow = 2)


