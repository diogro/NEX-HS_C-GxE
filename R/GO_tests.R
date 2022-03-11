source(here::here("R/go_functions.R"))

# en_head = makeEnrichment("data/output/SBM/clustering/head_weights-spearman_fdr-1e-02_mcmc_mode_hierarchical-SBM_gene-blocks")
# en_body = makeEnrichment("data/output/SBM/clustering/body_weights-spearman_fdr-1e-03_mcmc_mode_hierarchical-SBM_gene-blocks")
# en_head$WGCNA = llply(0:8, getEnrichment,
#                       folder_path="data/output/SBM/clustering/head_weights-spearman_fdr-1e-04_mcmc_mode_hierarchical-SBM_gene-blocks",
#                       clustering = "WGCNA")
# en_body$WGCNA = llply(0:9, getEnrichment,
#                       folder_path="data/output/SBM/clustering/body_weights-spearman_fdr-1e-05_mcmc_mode_hierarchical-SBM_gene-blocks",
#                       clustering = "WGCNA")
# saveRDS(en_head, here::here('data/enGo_head_fdr-1e-02.Rds'))
# saveRDS(en_body, here::here('data/enGo_body_fdr-1e-03.Rds'))

# en_head_abs = makeEnrichment("data/output/SBM/clustering/head_weights-spearman_fdr-1e-04_absolute_mcmc_mode_hierarchical-SBM_gene-blocks")
# en_body_abs = makeEnrichment("data/output/SBM/clustering/body_weights-spearman_fdr-1e-05_absolute_mcmc_mode_hierarchical-SBM_gene-blocks")
# saveRDS(en_head_abs, 'data/enGo_head_abs.Rds')
# saveRDS(en_body_abs, 'data/enGo_body_abs.Rds')

en_head = readRDS(here::here('data/enGo_head_fdr-1e-02.Rds'))
en_body = readRDS(here::here('data/enGo_body_fdr-1e-03.Rds'))

en_head_table = ldply(en_head$CP, function(x) x@result) %>%
  filter(p.adjust < 0.05, Count >= 4) %>% rename(Name = .id)
write.csv(en_head_table, file = here::here("data/output/SBM/GO/GOenrichment_head_fdr-1e-02.csv"))
en_body_table = ldply(en_body$CP, function(x) x@result) %>%
  filter(p.adjust < 0.05, Count >= 4) %>% rename(Name = .id)
write.csv(en_body_table, file = here::here("data/output/SBM/GO/GOenrichment_body_fdr-1e-03.csv"))

write.csv(en_head$summary, file = here::here("data/output/SBM/GO/Summary_head_fdr-1e-02.csv"))
write.csv(en_body$summary, file = here::here("data/output/SBM/GO/Summary_body_fdr-1e-03.csv"))



table_en = table(en_head$summary$n_enrich[en_head$summary$Nested_Level==1]!=0)
table_en
table_en/sum(table_en)

table_en = table(en_head$summary$n_enrich[en_head$summary$Nested_Level==2]!=0)
table_en
table_en/sum(table_en)

table_en = table(en_body$summary$n_enrich[en_body$summary$Nested_Level==1]!=0)
table_en
table_en/sum(table_en)

table_en = table(en_body$summary$n_enrich[en_body$summary$Nested_Level==2]!=0)
table_en
table_en/sum(table_en)



llply(en_body$summary$Name[en_body$summary$Nested_Level==2], CP_plot, en_body$CP, en_body$summary)

Level = 4
for(x in en_head$summary$Name[en_head$summary$Nested_Level==Level]){
  df = CP_print(x, enGo = en_head$CP, en_head$summary, recursive = FALSE, fdr = 1e-5)
  print(df)
}


wgcna_df_body =ldply(1:length(en_body$WGCNA), function(id){
  x <- en_body$WGCNA[[id]];
  dplyr::select(x@result, -geneID) %>%
    filter(p.adjust < 0.05, Count >= 4) %>%
    mutate(Module = id-1) %>%
    select(Module, everything())
})

wgcna_df_head =ldply(1:length(en_head$WGCNA), function(id){
      x <- en_head$WGCNA[[id]];
      dplyr::select(x@result, -geneID) %>%
        filter(p.adjust < 0.05, Count >= 4) %>%
        mutate(Module = id-1) %>%
        select(Module, everything())
})
wgcna_df_body %>% filter(grepl("cytoplasmic translation", Description))
wgcna_df_head %>% filter(grepl("ATP", Description))
cnetplot(en_head$WGCNA[[1]], node_label="category",showCategory = 10)
cnetplot(en_body$WGCNA[[5]], node_label="category",showCategory = 10)

wgcna_df %>% filter(Module==4)

translate_head = inner_join(en_head_table  %>% filter(grepl("cytoplasmic translation", Description)) %>%
                              select(-geneID, -pvalue),
           en_head$summary, by="Name") %>% filter(Nested_Level == 2)
neuro_head = inner_join(en_head_table  %>% filter(grepl("G protein-coupled receptor", Description)) %>%
                              select(-geneID, -pvalue),
                            en_head$summary, by="Name") %>% filter(Nested_Level == 2)

inner_join(en_head_table %>%
             select(-geneID, -pvalue),
           en_head$summary, by="Name") %>% filter(Nested_Level == 2 & Block %in% c(5))

inner_join(en_body_table %>%
             select(-geneID, -pvalue),
           en_body$summary, by="Name") %>% filter(Nested_Level == 2 & Block %in% c(12))

translate_body = inner_join(en_body_table  %>% filter(grepl("cytoplasmic translation", Description)) %>%
                              select(-geneID, -pvalue),
           en_body$summary, by="Name") %>% filter(Nested_Level == 2)
respiration_body = inner_join(en_body_table  %>% filter(grepl("fatty acid", Description)) %>%
                              select(-geneID, -pvalue),
                            en_body$summary, by="Name") %>% filter(Nested_Level == 2)

mean(c(translate_head$Assortatitvity, translate_body$Assortatitvity))
min(c(translate_head$Assortatitvity, translate_body$Assortatitvity))
max(c(translate_head$Assortatitvity, translate_body$Assortatitvity))


names_head = translate_head$Name
names_body = translate_body$Name
SBM_translation_list = llply(en_head$CP[names_head], cnetplot, node_label="category",showCategory = 10)
SBM_translation_list = append(SBM_translation_list, llply(en_head$CP[names_body], cnetplot, node_label="category",showCategory = 10))
plot_translation= plot_grid(plotlist = SBM_translation_list, ncol = 2,
                            labels = c(paste("Head -", names_head), paste("Body -", names_head)), scale = 0.9)
save_plot("~/Dropbox/labbio/articles/NEX_BodyHead_Control-SBM/figures/Translation_go_map.png", plot_translation, base_height = 5,
          base_asp = 2, ncol = 2, nrow = 4)

  inner_join(en_body_table  %>% filter(grepl("11-1-1", .id)) %>% select(-geneID, -pvalue) %>% rename(Name=.id),
             en_head$summary, by="Name") %>% filter(Nested_Level == 1)
getChild("11-1-1", en_body$summary)
cnetplot(en_body$CP$`11-1-1`, node_label="category",showCategory = 10)
translation_names = getChild("11-1-1", en_body$summary)[-1]
plot_11_1_1 = plot_grid(plotlist = llply(en_body$CP[translation_names], cnetplot,
                                         node_label="category",showCategory = 10),
                        labels = translation_names, ncol = 2, scale = 0.9)
save_plot("~/Dropbox/labbio/articles/NEX_BodyHead_Control-SBM/figures/plot_11_1_1.png", plot_11_1_1, base_height = 5,
          base_asp = 1.5, ncol = 2, nrow = 3)
translation_names = getChild("1-1-1", en_body$summary)[-1]
plot_grid(plotlist = llply(en_body$CP[translation_names], cnetplot, node_label="category",showCategory = 10))

names = getChild("0-0-0", en_head$summary)[-1]
SBM_neuro_list = llply(en_head$CP[names], cnetplot, node_label="category",showCategory = 10)
WGCNA_neuro = cnetplot(en_head$WGCNA[[5]], node_label="category",showCategory = 10)
WGCNA_photo = cnetplot(en_head$WGCNA[[7]], node_label="category",showCategory = 10)

neuro_list = SBM_neuro_list; neuro_list[[length(neuro_list)+1]] = WGCNA_neuro
plot_000 = plot_grid(plotlist = neuro_list, ncol = 2, labels = c(names, "WGCNA module 4"), scale = 0.9)
save_plot("~/Dropbox/labbio/articles/NEX_BodyHead_Control-SBM/figures/000_go_map.png", plot_000, base_height = 5,
          base_asp = 1.5, ncol = 2, nrow = 3)
plot_grid(plotlist = llply(en_head$CP[getChild("12-0-0", en_head$summary)[c(-1, -3, -5)]], cnetplot))

getChild("1-1", en_head$summary)[c(-1)]
plot_grid(plotlist = llply(en_head$CP_simple[getChild("1-1", en_head$summary)[c(-1)]], emplot))


cnetplot(en_head$CP$`9-4-0`)
cnetplot(en_head$CP$`1-1-1-1`)

cnetplot(en_head$CP$`30-2-2-1`)

cnetplot(en_head$WGCNA[[5]], node_label="category", showCategory = 20)
dotplot(en_head$WGCNA[[5]])

?dotplot


rbind(read_csv("data/output/SBM/clustering/head_weights-spearman_fdr-1e-04_mcmc_mode_hierarchical-SBM.csv") %>%
  mutate(tissue = "head"),
read_csv("data/output/SBM/clustering/body_weights-spearman_fdr-1e-05_mcmc_mode_hierarchical-SBM.csv") %>%
  mutate(tissue = "body")) %>% select(Gene, Degree, E_corr, B1:B4, tissue) %>%
  write_csv("block_partition-HeadBody-Control.csv")


bp2 <- clusterProfiler::simplify(en_head$WGCNA[[5]])
emplot(clusterProfiler::simplify(en_head$WGCNA[[2]]))
emplot(clusterProfiler::simplify(en_head$CP$`1-1-1`))

emplot(clusterProfiler::simplify(en_head$WGCNA[[3]]))
emplot(clusterProfiler::simplify(en_head$WGCNA[[4]]))
emplot(clusterProfiler::simplify(en_head$WGCNA[[5]]))
emplot(clusterProfiler::simplify(en_head$WGCNA[[6]]))
emplot(clusterProfiler::simplify(en_head$WGCNA[[7]]))
emplot(clusterProfiler::simplify(en_head$WGCNA[[8]]))

enrichplot::treeplot(enrichplot::pairwise_termsim(en_head$WGCNA[[5]]))
enrichplot::treeplot(enrichplot::pairwise_termsim(en_head$CP$`3-0`))

enrichplot::treeplot(enrichplot::pairwise_termsim(en_head$CP$`1-1-1`))
enrichplot::treeplot(enrichplot::pairwise_termsim(en_head$CP$`2-1`))
enrichplot::treeplot(enrichplot::pairwise_termsim(en_head$WGCNA[[2]]))



plot_grid(
emapplot(enrichplot::pairwise_termsim(en_head$WGCNA[[5]])),
empplot(enrichplot::pairwise_termsim(en_head$CP$`0-0-0-0`)))

emplot = function(x) emapplot(enrichplot::pairwise_termsim(x))

pak::pkg_install("GuangchuangYu/enrichplot")

auto_barplot = function(x, en = en_head$CP){
  x = en[[x]]
  if(nrow(summary(x)) == 0) return(NULL)
  current_go = clusterProfiler::simplify(x)
  pairwise_go <- pairwise_termsim(current_go) 
  barplot(pairwise_go)
}

level1 = getChild("6-0-0-0", en_head$summary)
pl = llply(level1, auto_barplot, en_head$CP)
plot_grid(plotlist = pl)

level1 = getChild("3-1-2-1", en_head$summary)
pl = llply(level1, auto_barplot, en_head$CP)
plot_grid(plotlist = pl)

level1 = getChild("8-4-1-1", en_head$summary)
pl = llply(level1, auto_barplot, en_head$CP)
plot_grid(plotlist = pl)
