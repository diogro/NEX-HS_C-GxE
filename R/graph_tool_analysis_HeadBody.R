if(!require(plyr)){install.packages("plyr"); library(plyr)}
if(!require(tidyverse)){install.packages("tidyverse"); library(tidyverse)}
if(!require(cowplot)){install.packages("cowplot"); library(cowplot)}
if(!require(ggthemes)){install.packages("ggthemes"); library(ggthemes)}
if(!require(xkcd)){install.packages("xkcd"); library(xkcd)}
if(!require(extrafont)){install.packages("extrafont"); library(extrafont)}
if(!require(superheat)){install.packages("superheat"); library(superheat)}
if(!require(mcclust)){install.packages("mcclust"); library(mcclust)}
if(!require(patchwork)){install.packages("patchwork"); library(patchwork)}

head_files = dir("data/output/SBM/clustering/", pattern = "head_cutoff-spearman_val-", full.names = T)
head_labels_raw = dir("data/output/SBM/clustering/", pattern = "head_cutoff-spearman_val-")

body_files = dir("data/output/SBM/clustering/", pattern = "body_cutoff-spearman_val-", full.names = T)
body_labels_raw = dir("data/output/SBM/clustering/", pattern = "body_cutoff-spearman_val-")


readClusteringFolder <- function(tissue_files, tissue_labels_raw) {
  tissue_labels_df = data.frame(cut0ff = sapply(str_split(tissue_labels_raw,"_"), function(x) str_split(x[[3]], "-")[[1]][2]),
             type = sapply(str_split(tissue_labels_raw,"_"), function(x) str_split(x[[4]], "-")[[1]][1])) %>%
    mutate(type =  ifelse(type == "non", "non.nested", "nested"),
           label = paste(type, cut0ff, sep = "."))
  
  readClustering = function(file){
    x = read_csv(file) %>%
      mutate_at(vars(matches("B")), ~as.numeric(factor(.)))
  }
  
  tissue_labels_df$cut0ff
  tissue_clusters_raw = lapply(tissue_files, readClustering)
  tissue_clusters = list(nested = tissue_clusters_raw[tissue_labels_df$type == "nested"],
                         non.nested =  tissue_clusters_raw[tissue_labels_df$type == "non.nested"])
  names(tissue_clusters[['nested']]) = tissue_labels_df$cut0ff[tissue_labels_df$type == "nested"]
  names(tissue_clusters[['non.nested']]) = tissue_labels_df$cut0ff[tissue_labels_df$type == "non.nested"]
  
  tissue_clusters
}

head_clusters = readClusteringFolder(head_files, head_labels_raw)
body_clusters = readClusteringFolder(body_files, body_labels_raw)

n_clusters_df = rbind(
  data.frame(n_clusters = sapply(head_clusters[['non.nested']], function(x) length(unique((x$Block)))),
             n_genes = sapply(head_clusters[['non.nested']], function(x) length(unique((x$Gene)))),
             cutOff = as.numeric(names(head_clusters[['non.nested']])),
             type = "non.nested", tissue = "head"),
  data.frame(n_clusters = sapply(head_clusters[['nested']], function(x) length(unique((x$B1)))),
             n_genes = sapply(head_clusters[['nested']], function(x) length(unique((x$Gene)))),
             cutOff = as.numeric(names(head_clusters[['nested']])),
             type = "nested", tissue = "head") ,
  data.frame(n_clusters = sapply(body_clusters[['non.nested']], function(x) length(unique((x$Block)))),
             n_genes = sapply(body_clusters[['non.nested']], function(x) length(unique((x$Gene)))),
             cutOff = as.numeric(names(body_clusters[['non.nested']])),
             type = "non.nested", tissue = "body"  ),
  data.frame(n_clusters = sapply(body_clusters[['nested']], function(x) length(unique((x$B1)))),
             n_genes = sapply(body_clusters[['nested']], function(x) length(unique((x$Gene)))),
             cutOff = as.numeric(names(body_clusters[['nested']])),
             type = "nested", tissue = "body")  
                )

n_clusters_plot = plot_grid(n_clusters_df %>%
  filter(tissue == "head") %>%
  ggplot(aes(cutOff, n_clusters, groups = type, color= type)) + 
    geom_line(size = 1) + 
    theme_tufte() + background_grid() + theme(legend.position = "top") +
    ggtitle("Number of clusters - Head - Nested vs non-nested SBM"),
n_clusters_df %>%
  filter(tissue == "body") %>%
  ggplot(aes(cutOff, n_clusters, groups = type, color= type)) + 
  geom_line(size = 1) + 
  theme_tufte() + background_grid() + theme(legend.position = "top") +
  ggtitle("Number of clusters - Body - Nested vs non-nested SBM"), ncol = 2)
save_plot(filename = "data/output/SBM/plots/number_of_gene_clusters_nested-vs-non-nested.png", 
          n_clusters_plot, base_height = 5, ncol = 2, base_asp = 1.2)

levels_array = laply(head_clusters[['nested']], function(x) 
  laply(select(x, matches("B")), 
        function(df) length(unique(df))))    
levels_array_df = data.frame(levels_array, cutOff = as.numeric(names(head_clusters[['nested']])))
names(levels_array_df) =  gsub("X", "Level.", names(levels_array_df))
nested_plot_head = levels_array_df %>% 
  pivot_longer(Level.1:Level.5) %>%
  ggplot(aes(cutOff, value, group = name, color = name)) +
  geom_line(size = 1) + 
  theme_tufte() + background_grid() + theme(legend.position = "top") +
  ggtitle("Number of clusters - Head - Nested - All levels")


levels_array = laply(body_clusters[['nested']], function(x) 
  laply(select(x, matches("B")), 
        function(df) length(unique(df))))    
levels_array_df = data.frame(levels_array, cutOff = as.numeric(names(body_clusters[['nested']])))
names(levels_array_df) =  gsub("X", "Level.", names(levels_array_df))
nested_plot_body = levels_array_df %>% 
  pivot_longer(Level.1:Level.5) %>%
  ggplot(aes(cutOff, value, group = name, color = name)) +
  geom_line(size = 1) + 
  theme_tufte() + background_grid() + theme(legend.position = "top") +
  ggtitle("Number of clusters - Body - Nested - All levels")

nested_plot = plot_grid(nested_plot_head, nested_plot_body, ncol = 2)
save_plot(filename = "data/output/SBM/plots/number_of_gene_clusters_nested-all_levels.png", 
          nested_plot, base_height = 5, ncol = 2, base_asp = 1.2)
