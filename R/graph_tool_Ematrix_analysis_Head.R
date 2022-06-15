if(!require(plyr)){install.packages("plyr"); library(plyr)}
if(!require(tidyverse)){install.packages("tidyverse"); library(tidyverse)}
if(!require(cowplot)){install.packages("cowplot"); library(cowplot)}
if(!require(ggthemes)){install.packages("ggthemes"); library(ggthemes)}
if(!require(xkcd)){install.packages("xkcd"); library(xkcd)}
if(!require(extrafont)){install.packages("extrafont"); library(extrafont)}
if(!require(superheat)){install.packages("superheat"); library(superheat)}
if(!require(mcclust)){install.packages("mcclust"); library(mcclust)}
if(!require(patchwork)){install.packages("patchwork"); library(patchwork)}

header = "head_weights-spearman_fdr-1e-02_mcmc_mode"

getBlockSizedf = function(level, block_df, all = FALSE, draw_level=level){
  upper = paste0("B", level+1)
  if(all){
    s_block_df = unique(block_df[,paste0("B", draw_level:(level+1))])
  } else
    s_block_df = unique(block_df[,paste0("B", level:(level+1))])

  t_upper = table(s_block_df[,upper])
  block_sizes = t_upper[match(unique(as.data.frame(s_block_df)[,upper]), names(t_upper))]
  sum_block_sizes = cumsum(block_sizes)
  b_size_mat = matrix(NA, nrow(block_sizes), 2)
  for(i in 1:nrow(block_sizes)){
    if(i == 1)
      b_size_mat[i,] = c(0.5, sum_block_sizes[i]+0.5)
    else
      b_size_mat[i,] = c(sum_block_sizes[i-1]+0.5, sum_block_sizes[i]+0.5)
  }
  b_size_df = data.frame(start = b_size_mat[,1],
                         end = b_size_mat[,2])
  return(b_size_df)
}
makeEmatrixPlots = function(header,
                            data_path = "data/output/SBM/clustering/",
                            plot_path = "data/output/SBM/plots/",
                            levels = 5){
  block_df = read_csv(file.path(data_path, paste0(header, "_hierarchical-SBM.csv")))
  block_summary = read_csv(file.path(data_path, paste0(header, "_hierarchical-SBM_gene-blocks/block_summary.csv"))) %>%
    filter(Nested_Level == 1) %>% select(Block, Internal_degree, Assortatitvity) %>%
    rename(B1 = Block)
  block_df = inner_join(block_df, block_summary)
  block_df = block_df %>%
    arrange(B6, B5, B4, B3, B2, desc(Assortatitvity)) %>%
    select(-1)

  e_mats = vector("list", levels)
  for (i in seq_along(e_mats)){
    e_matrix = read_csv(file.path(data_path, paste0(header, "_hierarchical-SBM_e_matrix_level",i-1,".csv")))
    e_matrix[,1] = NULL
    e_matrix = as.matrix(e_matrix)
    rownames(e_matrix) = colnames(e_matrix) = as.character(colnames(e_matrix))
    blocks = as.character(unique(as.data.frame(block_df)[,paste0("B", i)]))
    e_matrix= e_matrix[blocks, blocks]
    #e_matrix = e_matrix[-61, -61]
    diag(e_matrix) = diag(e_matrix)/2
    e_mats[[i]] = e_matrix
  }
  makePlot = function(level, e_mats, block_df){
    e_matrix = e_mats[[level]]
    colnames(e_matrix) = rownames(e_matrix) = paste0("b", rownames(e_matrix))
    melted_cormat <- reshape2::melt(e_matrix)
    melted_cormat[melted_cormat == 0] = NA
    #plot heatmap
    plot = ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) +
      geom_tile() +
      scale_fill_viridis_c(alpha = 1)  + theme_cowplot() +
      labs(y = "Blocks", x = "Blocks") + ggtitle(paste("Level", level)) +
      theme(axis.text.x = element_text(angle=45, vjust=0.6), legend.position = "none")
    if(level < levels){
      for(i in level:(levels-1)){
        b_size_df = getBlockSizedf(i, block_df, all = TRUE, level)
        plot = plot + geom_rect(data = b_size_df, color = "tomato3", alpha = 0, size = 1,
                                aes(x = NULL, y = NULL, fill = NULL, xmin=start, xmax=end,
                                    ymin=start, ymax=end))
      }
    }
    return(plot)
  }
  plot_list = lapply(seq_along(1:levels), makePlot, e_mats, block_df)

  all_plots = plot_grid(plotlist = plot_list)
  save_plot(file.path(plot_path, paste0(header, "_E_matrices.png")), all_plots,
            base_height = 10, base_asp = 1.2, ncol = 3, nrow = 2)
  return(list(df = block_df, E = e_mats, plots = all_plots, plot_list = plot_list))
}

out_fdr_1e2_head = makeEmatrixPlots("head_weights-spearman_fdr-1e-02_mcmc_mode", levels = 4)
out_fdr_1e3_body = makeEmatrixPlots("body_weights-spearman_fdr-1e-03_mcmc_mode", levels = 4)

library(patchwork)
pak::pkg_install("ggeasy")
library(ggeasy)


img1 <- magick::image_read("data/output/SBM/guide_plots/body/fdr-1e-03/full_plot.png")
body_full = ggplot2::ggplot() + ggplot2::annotation_custom(grid::rasterGrob(img1,
                                               width=ggplot2::unit(1,"npc"),
                                               height=ggplot2::unit(1,"npc")),
                               -Inf, Inf, -Inf, Inf) + ggtitle("C. Body - Full nested SBM") + 
                               theme_cowplot() + easy_remove_axes()
img2 <- magick::image_read("data/output/SBM/guide_plots/head/fdr-1e-02/full_plot.png")
head_full = ggplot2::ggplot() + ggplot2::annotation_custom(grid::rasterGrob(img2,
                                               width=ggplot2::unit(1,"npc"),
                                               height=ggplot2::unit(1,"npc")),
                               -Inf, Inf, -Inf, Inf) + ggtitle("D. Head - Full nested SBM") + 
                               theme_cowplot() + easy_remove_axes()


plot = out_fdr_1e3_body$plot_list[[1]] + ggtitle("A. Body - SBM Level-1") + 
       out_fdr_1e2_head$plot_list[[1]] + ggtitle("B. Head - SBM Level-1") +
        body_full + head_full
plot
save_plot("~/Dropbox/labbio/articles/NEX_BodyHead_Control-SBM/figures/SBM_Ematrix.png", plot, base_height = 10, ncol = 2, nrow = 2, base_asp = 1.)

