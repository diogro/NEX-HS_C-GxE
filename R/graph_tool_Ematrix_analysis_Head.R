if(!require(plyr)){pak::pkg_install("plyr"); library(plyr)}
if(!require(tidyverse)){pak::pkg_install("tidyverse"); library(tidyverse)}
if(!require(cowplot)){pak::pkg_install("cowplot"); library(cowplot)}
if(!require(ggthemes)){pak::pkg_install("ggthemes"); library(ggthemes)}
if(!require(extrafont)){pak::pkg_install("extrafont"); library(extrafont)}
if(!require(superheat)){pak::pkg_install("superheat"); library(superheat)}
if(!require(mcclust)){pak::pkg_install("mcclust"); library(mcclust)}
if(!require(patchwork)){pak::pkg_install("patchwork"); library(patchwork)}

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
    filter(Nested_Level == 1) %>% select(Block, Internal_degree, Assortativity) %>%
    rename(B1 = Block)
  block_df = inner_join(block_df, block_summary)
  block_df = block_df %>%
    arrange(B6, B5, B4, B3, B2, desc(Assortativity)) %>%
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
  makePlotE = function(level, e_mats, block_df){
    e_matrix = e_mats[[level]]
    colnames(e_matrix) = rownames(e_matrix) = paste0("b", rownames(e_matrix))
    melted_cormat <- reshape2::melt(e_matrix)
    melted_cormat[melted_cormat == 0] = NA
    #plot heatmap
    plot = ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) +
      geom_tile() +
      scale_fill_viridis_c(alpha = 1)  + 
      theme_cowplot() +
      labs(y = "Level-1 Blocks", x = "Level-1 Blocks") + 
      ggtitle(paste("Level", level)) +
      scale_y_discrete(label = gsub("b", "", colnames(e_matrix))) + 
      scale_x_discrete(label = gsub("b", "", colnames(e_matrix)), guide = guide_axis(n.dodge = 2)) + 
      theme(axis.text.x = element_text(angle=35, vjust=0.3), legend.position = "bottom", 
                                      legend.key.width= unit(22, 'mm'), legend.key.height= unit(2.2, 'mm'), legend.title = element_blank())
    if(level < levels){
      for(i in level:(levels-1)){
        b_size_df = getBlockSizedf(i, block_df, all = TRUE, level)
        plot = plot + geom_rect(data = b_size_df, color = "tomato3", alpha = 0, size = 0.5,
                                aes(x = NULL, y = NULL, fill = NULL, xmin=start, xmax=end,
                                    ymin=start, ymax=end))
      }
    }
    return(plot)
  }
  makePlotDeg = function(level, block_df){
    block_df[[paste0("B", level)]] = factor(block_df[[paste0("B", level)]], 
                                            level = unique(block_df[[paste0("B", level)]]))
    plot = ggplot(data = block_df, aes(x=B1, y=Degree)) + geom_boxplot() + 
        scale_y_continuous(breaks = c(0, 1000, 2000), minor_breaks = c(500, 1500)) +
        coord_flip() + theme_cowplot() + background_grid(minor = "xy") + easy_remove_axes() + 
        theme(axis.title.x = element_text(), axis.text.x = element_text())
    return(plot)
  }
  plot_list = lapply(seq_along(1:levels), makePlotE, e_mats, block_df)
  plot_list_deg = lapply(seq_along(1:levels), makePlotDeg, block_df)
  all_plots = plot_grid(plotlist = plot_list)
  save_plot(file.path(plot_path, paste0(header, "_E_matrices.png")), all_plots,
            base_height = 10, base_asp = 1.2, ncol = 3, nrow = 2)
  return(list(df = block_df, E = e_mats, plots = all_plots, plot_list = plot_list, plot_list_deg = plot_list_deg))
}

out_fdr_1e2_head = makeEmatrixPlots("head_weights-spearman_fdr-1e-02_mcmc_mode", levels = 5)
out_fdr_1e3_body = makeEmatrixPlots("body_weights-spearman_fdr-1e-03_mcmc_mode", levels = 5)

library(patchwork)
pak::pkg_install("ggeasy")
library(ggeasy)

{
  title_size = 10
  axis_size = 8
  tick_size = 4
  img1 <- magick::image_read("data/output/SBM/guide_plots/body/fdr-1e-03/trial.png")
  body_full = ggplot2::ggplot() + ggplot2::annotation_custom(grid::rasterGrob(img1,
                                                width=ggplot2::unit(1,"npc"),
                                                height=ggplot2::unit(0.975,"npc")),
                                -Inf, Inf, -Inf, Inf) + ggtitle("C.") + 
                                theme_cowplot() + easy_remove_axes() + theme(plot.title = element_text(size = title_size))
  img2 <- magick::image_read("data/output/SBM/guide_plots/head/fdr-1e-02/trial.png")
  head_full = ggplot2::ggplot() + ggplot2::annotation_custom(grid::rasterGrob(img2,
                                                width=ggplot2::unit(1,"npc"),
                                                height=ggplot2::unit(1,"npc")),
                                -Inf, Inf, -Inf, Inf) + ggtitle("D.") + 
                                theme_cowplot() + easy_remove_axes() + theme(plot.title = element_text(size = title_size))

  plot = out_fdr_1e3_body$plot_list[[1]] + ggtitle("A. Body") + theme(plot.title = element_text(size = title_size), 
                                                                      axis.title.x = element_text(size = axis_size),
                                                                      axis.title.y = element_text(size = axis_size),
                                                                      axis.text.y = element_blank(),
                                                                      axis.text.x = element_text(size = tick_size, 
                                                                                    margin = margin(t = 0, r = 0, b = 0.2, l = 0)),
                                                                      axis.ticks.x = element_line(size = .2),
                                                                      axis.ticks.length=unit(.05, "cm"),
                                                                      axis.ticks.y = element_blank(),
                                                                      axis.line = element_blank(),
                                                                      legend.key.width= unit(1.5, 'cm'),
                                                                      legend.text = element_text(size = tick_size+1)) + 
        out_fdr_1e2_head$plot_list[[1]] + ggtitle("B. Head") + theme(plot.title = element_text(size = title_size), 
                                                                      axis.title.x = element_text(size = axis_size),
                                                                      axis.title.y = element_text(size = axis_size),
                                                                      axis.text.y = element_blank(),
                                                                      axis.text.x = element_text(size = tick_size, 
                                                                                    margin = margin(t = 0, r = 0, b = 0.2, l = 0)),
                                                                      axis.ticks.x = element_line(size = .2),
                                                                      axis.ticks.length=unit(.05, "cm"),
                                                                      axis.ticks.y = element_blank(),
                                                                      axis.line = element_blank(),
                                                                      legend.key.width= unit(1.5, 'cm'),
                                                                      legend.text = element_text(size = tick_size+1)) + 
          body_full + head_full
  save_plot("~/Dropbox/labbio/articles/SBM_manuscript/figures/SBM_Ematrix.png", plot, 
            base_width = 7.5, base_height = 7.5)
}


{layout <- "AAAAABCCCCCD
AAAAABCCCCCD
AAAAABCCCCCD
AAAAABCCCCCD
AAAAABCCCCCD
"}
plot = out_fdr_1e3_body$plot_list[[1]] + ggtitle("A. Body")  + 
            out_fdr_1e3_body$plot_list_deg[[1]] + plot_layout(design = layout) +
            out_fdr_1e2_head$plot_list[[1]] + ggtitle("B. Head")  + 
            out_fdr_1e2_head$plot_list_deg[[1]] + plot_layout(design = layout)
plot
save_plot("~/Dropbox/labbio/articles/NEX_BodyHead_Control-SBM/figures/SBM_Ematrix_with_degree.png", plot, 
          base_height = 11, ncol = 2, nrow = 1, base_asp = 1. + 1/6)
