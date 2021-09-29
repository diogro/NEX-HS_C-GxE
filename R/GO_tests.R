# BiocManager::install("org.Dm.eg.db")
# BiocManager::install("hfang-bristol/XGR", dependencies=T)

library(XGR)
library(plyr)
library(cowplot)
library(ggrepel)
library("AnnotationDbi")
library("org.Dm.eg.db")
# flyGO <- xRDataLoader(RData.customised = 'org.Dm.egGOBP',
#                       RData.location = "https://github.com/hfang-bristol/RDataCentre/blob/master/dnet/1.0.7")
# saveRDS(flyGO, "org.Dm.egGOBP.2021-09-27")
flyGO = readRDS("org.Dm.egGOBP.2021-09-27")

#BiocManager::install("clusterProfiler")
library(clusterProfiler)

getLevel = function(x, n_levels = 4){
  x = gsub(".csv", "", x)
  x = gsub("/", "", x)
  n_levels - length(strsplit(x, "-")[[1]]) + 1
}
getChild = function(parent, b_df = block_summary){
  if(!parent %in% b_df$Name) stop("Parent block not found.")
  childs = b_df %>%
    filter(Nested_Level == getLevel(parent)-1, Parent == strsplit(parent, "-")[[1]][1]) %>%
    as.data.frame %>% `[[`("Name")
  c(parent, childs)
}

getGeneList = function(block_number, level, folder_path, file_name = NULL){
  if(is.null(file_name)){
    file_path = file.path(folder_path, paste0("Level_", level))
    file = dir(file_path, pattern = paste0("^", block_number), full.names = T)
  } else{
    level = getLevel(file_name)
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
                         folder_path, file_name = NULL, type = c("clusterProfiler", "XGR")){
  type = match.arg(type)
  x = getGeneList(block_number, level, folder_path, file_name)
  if(type == "clusterProfiler")
  enGo = enrichGO(gene          = x$en$genes,
                  universe      = x$en$background,
                  OrgDb         = org.Dm.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  else
    enGo = xEnricherGenes(data = x$en$genes,
                         ontology.customised = flyGO,
                         ontology = 'NA',
                         background = x$en$background,
                         verbose = TRUE,
                         check.symbol.identity = TRUE,
                         ontology.algorithm = "none")
  enGo
}
XGR_plot = function(x){
  l = enGo_XGR[getChild(x)]
  l = l[!sapply(l, is.null)]
  if(length(l) > 2){
    print(x)
    goplot <- xEnrichCompare(l, displayBy="fc", FDR.cutoff = 0.05, wrap.width = 45) +
      scale_fill_brewer(palette='Set2') +
      ggtitle('GO Enrichment Biological Process (FDR < 0.05)')
  } else
    goplot = NULL
  return(goplot)
}
CP_plot = function(x){
  l = enGo_CP_simple[getChild(x)]
  l = l[!sapply(l, is.null)]
  if(length(l) > 2){
    print(x)
    df = ldply(l, function(x) select(x@result, -geneID), .id = "Block") %>%
      filter(p.adjust < 0.05)
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
CP_print = function(x, fdr = 0.05, simple = FALSE, recursive=TRUE){
  if(simple){
    enGo = enGo_CP_simple
  } else{
    enGo = enGo_CP
  }
  if(recursive){
    out = ldply(enGo[getChild(x)],
          function(d) d@result %>%
            dplyr::select(-geneID) %>%
            filter(p.adjust < fdr), .id =  "Block")
  } else
    out = enGo[[x]]@result %>%
      dplyr::select(-geneID) %>%
      filter(p.adjust < fdr) %>%
      mutate(Block = x) %>%
      select(Block, everything())
  return(out)
}

block_path = "data/output/SBM/clustering/head_weights-spearman_fdr-1e-06_mcmc_mode_hierarchical-SBM_gene-blocks"
level = 1
block  = 11
block_summary = read.csv(file.path(block_path, "block_summary.csv"))
block_summary$Name = block_summary$File %>% {gsub(".csv", "", .)} %>% {gsub("/", "", .)}
block_summary$Parent = block_summary$Name %>%
  llply(strsplit, split="-") %>% sapply(`[[`, 1) %>% sapply(`[`, 2)

enGo_XGR = llply(block_summary$File,
              function(file)
                getEnrichment(file_name = file,
                              folder_path = block_path,
                              type = "XGR"), .progress = "text")
names(enGo_XGR) = block_summary$Name
enGo_CP = llply(block_summary$File,
                  function(file)
                    getEnrichment(file_name = file,
                                  folder_path = block_path,
                                  type = "clusterProfiler"), .progress = "text")
names(enGo_CP) = block_summary$Name

enGo_CP_simple = lapply(enGo_CP, clusterProfiler::simplify, cutoff=0.7, by="p.adjust", select_fun=min)

x = enGo_XGR[[2]]
x
dplyr::select(enGo_CP[[2]]@result, -geneID) %>% head()

block_summary$p.adjust = laply(enGo_CP, function(x) select(x@result, p.adjust)[1,])
block_summary$n_enrich = laply(enGo_CP,
                               function(x) select(x@result, p.adjust) %>%
                                 filter(p.adjust < 0.05) %>% nrow)
block_summary$n_enrich_simple = laply(enGo_CP_simple,
                               function(x) select(x@result, p.adjust) %>%
                                 filter(p.adjust < 0.05) %>% nrow)

table_en = table(block_summary$n_enrich[block_summary$Nested_Level==1]!=0)
table_en/sum(table_en)

ggplot(filter(block_summary, Nested_Level == 1), aes(Assortatitvity, -log2(p.adjust), color = Parent)) +
  geom_point() + geom_label_repel(aes(label = Block))
ggplot(filter(block_summary, Nested_Level == 1), aes(N_genes, -log2(p.adjust), color = Parent)) +
  geom_point() + geom_label_repel(aes(label = Block))
ggplot(filter(block_summary, Nested_Level == 1), aes(Assortatitvity, n_enrich, color = Parent)) +
  geom_point() + geom_label_repel(aes(label = Block))
ggplot(filter(block_summary, Nested_Level == 1), aes(N_genes, n_enrich, color = Parent)) +
  geom_point() + geom_label_repel(aes(label = Block))
ggplot(filter(block_summary, Nested_Level == 1), aes(Assortatitvity, n_enrich_simple, color = Parent)) +
  geom_point() + geom_label_repel(aes(label = Block))
ggplot(filter(block_summary, Nested_Level == 1), aes(N_genes, n_enrich_simple, color = Parent)) +
  geom_point() + geom_label_repel(aes(label = Block))

ggplot(filter(block_summary, Nested_Level == 2), aes(Assortatitvity, -log2(p.adjust), color = Parent)) +
  geom_point() + geom_label_repel(aes(label = Block))
ggplot(filter(block_summary, Nested_Level == 2), aes(N_genes, -log2(p.adjust), color = Parent)) +
  geom_point() + geom_label_repel(aes(label = Block))
ggplot(filter(block_summary, Nested_Level == 2), aes(Assortatitvity, n_enrich/N_genes, color = Parent)) +
  geom_point() + geom_label_repel(aes(label = Block))+ geom_smooth(method = lm, color = "black")
ggplot(filter(block_summary, Nested_Level == 2), aes(N_genes, n_enrich, color = Parent)) +
  geom_point() + geom_label_repel(aes(label = Block))
ggplot(filter(block_summary, Nested_Level == 2), aes(Assortatitvity, n_enrich_simple, color = Parent)) +
  geom_point() + geom_label_repel(aes(label = Block))
ggplot(filter(block_summary, Nested_Level == 2), aes(N_genes, n_enrich_simple, color = Parent)) +
  geom_point() + geom_label_repel(aes(label = Block))


goplot_list = llply(block_summary$Name[block_summary$Nested_Level==3], XGR_plot)
save_plot("go_head_level_11-1-0-super_translation.png", XGR_plot("11-1-0"), base_height = 7, base_asp = 0.25, ncol=4)

l = enGo_CP_simple[block_summary$Name[block_summary$Nested_Level==3]]

llply(block_summary$Name[block_summary$Nested_Level==2], CP_plot)

save_plot("go_head_level_4-0-0_gueto.png", goplot, base_height = 9, base_asp = 0.25, ncol=5)

Level = 3
for(x in block_summary$Name[block_summary$Nested_Level==Level]){
  df = CP_print(x)
  print(df)
}
