#pak::pkg_install(c("ggrepel", "AnnotationDbi", "org.Dm.eg.db", 
#                   "clusterProfiler", "ggnewscale", "enrichplot", 
#                   "tidyverse", "plyr",  "cowplot", "ggrepel", 
#                   "ggpubr", "ggthemes"))
library(plyr)
library(cowplot)
library(ggrepel)
library(tidyverse)

library(AnnotationDbi)
library(org.Dm.eg.db)
library(clusterProfiler)
library(ggnewscale)
library(enrichplot)
library(ggthemes)

getLevel = function(x, n_levels = 5){
  x = gsub(".csv", "", x)
  x = gsub("/", "", x)
  n_levels - length(strsplit(x, "-")[[1]]) + 1
}
getChild = function(parent, b_df){
  if(!parent %in% b_df$Name) stop("Parent block not found.")
  n_level = max(b_df$Nested_Level)
  childs = b_df %>%
    filter(Nested_Level == getLevel(parent, n_level)-1, Parent == strsplit(parent, "-")[[1]][1]) %>%
    as.data.frame %>% `[[`("Name")
  c(parent, childs)
}

getGeneList = function(block_number, level, folder_path, file_name = NULL, clustering){
  if(is.null(file_name)){
    if(clustering == "SBM"){
      file_path = file.path(folder_path, paste0("Level_", level))
      file = dir(file_path, pattern = paste0("^", block_number, "-"), full.names = T)
    } else if(clustering == "WGCNA"){
      file_path = file.path(folder_path, "WGCNA")
      file = dir(file_path, pattern = paste0("^wgcna_", block_number), full.names = T)
    }
  } else{
    block_summary = read.csv(file.path(folder_path, "block_summary.csv"))
    n_level = max(block_summary$Nested_Level)
    level = getLevel(file_name, n_level)
    file = file.path(folder_path, paste0("Level_", level), file_name)
  }
  gene_list = read.csv(file, header = FALSE)[,1]

  background = read.csv(file.path(folder_path, "background.csv"), header = FALSE)[,1]
  fb = list(genes = gene_list, background = background)
  en = list(genes = bitr(fb$genes, fromType="FLYBASE",
                         toType="ENTREZID", OrgDb="org.Dm.eg.db")$ENTREZID,
            background = bitr(fb$background, fromType="FLYBASE",
                              toType="ENTREZID", OrgDb="org.Dm.eg.db")$ENTREZID)
  list(fb = fb, en = en)
}
getEnrichment = function(block_number = NULL, level = NULL,
                         folder_path, file_name = NULL,
                         clustering = c("SBM", "WGCNA"),
                         type = c("clusterProfiler", "XGR", "group"),
                         ont = c("BP", "MF", "CC"),
                         go_level = 3){
  clustering = match.arg(clustering)
  type = match.arg(type)
  ont = match.arg(ont)
  x = getGeneList(block_number, level, folder_path, file_name, clustering = clustering)
  if(type == "clusterProfiler")
    enGo = enrichGO(gene          = x$en$genes,
                    universe      = x$en$background,
                    OrgDb         = org.Dm.eg.db,
                    ont           = ont,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
  else if(type == "group")
    enGo = groupGO(gene    = x$en$genes,
                   OrgDb   = org.Dm.eg.db,
                   ont = ont,
                   level = go_level,
                   readable = TRUE)
  else
    enGo = xEnricherGenes(data = x$en$genes,
                          ontology.customised = flyGO[[ont]],
                          ontology = NA,
                          background = x$en$background,
                          verbose = TRUE,
                          check.symbol.identity = TRUE,
                          ontology.algorithm = "none")
  enGo
}
XGR_plot = function(x, enGo, summary, ont = c("BP", "MF", "CC"), fdr = 0.05){
  ont = match.arg(ont)
  ontlogies = list(BP = 'Biological Process',
                   MF = 'Molecular Function',
                   CC = 'Cellular Component')
  title = paste0('GO Enrichment ', ontlogies[[ont]], ' (FDR < ', fdr, ')')
  l = enGo[getChild(x, summary)]
  l = l[!sapply(l, is.null)]
  if(length(l) >= 2){
    print(x)
    goplot <- xEnrichCompare(l, displayBy="fc", FDR.cutoff = fdr, wrap.width = 45) +
      scale_fill_brewer(palette='Set2') +
      ggtitle(title)
  } else
    goplot = NULL
  return(goplot)
}
CP_plot = function(x, enGo_CP_simple, summary){
  l = enGo_CP_simple[getChild(x, summary)]
  l = l[!sapply(l, is.null)]
  if(length(l) > 2){
    print(x)
    df = ldply(l, function(x) dplyr::select(x@result, -geneID), .id = "Block") %>%
      filter(p.adjust < 0.05, count >= 5)
    goplot <- df %>%
      ggplot(aes(Description, -log2(p.adjust))) +
      geom_bar(stat="identity")  + coord_flip()
    if(length(unique(df$Block))>1){
      goplot <- goplot + facet_wrap(~Block, nrow = 1)
    } else
      goplot <- ifelse(is.na(df$Block[1]), NULL, goplot + ggtitle(df$Block[1]))
  } else
    goplot = NULL
  return(goplot)
}
CP_print = function(x, enGo, summary, fdr = 0.05, recursive=TRUE){
  if(recursive){
    out = ldply(enGo[getChild(x, summary)],
                function(d) d@result %>%
                  dplyr::select(-geneID) %>%
                  filter(p.adjust < fdr), .id =  "Block")
  } else
    out = enGo[[x]]@result %>%
      dplyr::select(-geneID) %>%
      filter(p.adjust < fdr) %>%
      mutate(Block = x) %>%
      dplyr::select(Block, everything())
  return(out)
}

makeEnrichment = function(block_path){
  block_summary = read.csv(file.path(block_path, "block_summary.csv"))
  block_summary$Name = block_summary$File %>% {gsub(".csv", "", .)} %>% {gsub("/", "", .)}
  block_summary$Parent = block_summary$Name %>%
    llply(strsplit, split="-") %>% sapply(`[[`, 1) %>% sapply(`[`, 2)

  # enGo_XGR = llply(block_summary$File,
  #                  function(file)
  #                    getEnrichment(file_name = file,
  #                                  folder_path = block_path,
  #                                  type = "XGR"))
  # names(enGo_XGR) = block_summary$Name
  # enGo_XGR_CC = llply(block_summary$File,
  #                     function(file)
  #                       getEnrichment(file_name = file,
  #                                     folder_path = block_path,
  #                                     type = "XGR", ont = "CC"))
  # names(enGo_XGR_CC) = block_summary$Name
  enGo_CP = llply(block_summary$File,
                  function(file)
                    getEnrichment(file_name = file,
                                  folder_path = block_path,
                                  type = "clusterProfiler"))
  names(enGo_CP) = block_summary$Name

  enGo_CP_simple = lapply(enGo_CP, clusterProfiler::simplify, cutoff=0.7, by="p.adjust", select_fun=min)

  block_summary$p.adjust = laply(enGo_CP, function(x) dplyr::select(x@result, p.adjust)[1,])
  block_summary$n_enrich = laply(enGo_CP,
                                 function(x) dplyr::select(x@result, p.adjust, Count) %>%
                                   filter(p.adjust < 0.05, Count >= 4) %>% nrow)
  block_summary$n_enrich_simple = laply(enGo_CP_simple,
                                        function(x) dplyr::select(x@result, p.adjust, Count) %>%
                                          filter(p.adjust < 0.05, Count >= 4) %>% nrow)
  return(list(summary = block_summary,
              CP = enGo_CP,
              CP_simple = enGo_CP_simple
              # XGR = enGo_XGR,
              # XGR_CC = enGo_XGR_CC
              ))
}

go_plots = function(local_go_upper, cats = 20, ...){
  pw_upper <- pairwise_termsim(local_go_upper) 
  plot_upper <- emapplot(pw_upper, showCategory = cats, ...) + 
            theme_tufte() + 
            theme(legend.position = "none") + 
            theme(plot.title = element_text(size=28),
                  axis.title = element_blank(),
                  axis.text = element_blank(),
                  axis.ticks = element_blank()) 
  plot_upper$data$color = mean(plot_upper$data$color)
  plot_upper
}
