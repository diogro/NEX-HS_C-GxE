source("R/go_functions.R")

# en_head = makeEnrichment("data/output/SBM/clustering/head_weights-spearman_fdr-1e-04_mcmc_mode_hierarchical-SBM_gene-blocks")
# en_body = makeEnrichment("data/output/SBM/clustering/body_weights-spearman_fdr-1e-05_mcmc_mode_hierarchical-SBM_gene-blocks")
# en_head$WGCNA = llply(0:8, getEnrichment,
#                       folder_path="data/output/SBM/clustering/head_weights-spearman_fdr-1e-04_mcmc_mode_hierarchical-SBM_gene-blocks",
#                       clustering = "WGCNA")
# en_body$WGCNA = llply(0:9, getEnrichment,
#                       folder_path="data/output/SBM/clustering/body_weights-spearman_fdr-1e-05_mcmc_mode_hierarchical-SBM_gene-blocks",
#                       clustering = "WGCNA")
saveRDS(en_head, 'data/enGo_head.Rds')
saveRDS(en_body, 'data/enGo_body.Rds')

# en_head_abs = makeEnrichment("data/output/SBM/clustering/head_weights-spearman_fdr-1e-04_absolute_mcmc_mode_hierarchical-SBM_gene-blocks")
# en_body_abs = makeEnrichment("data/output/SBM/clustering/body_weights-spearman_fdr-1e-05_absolute_mcmc_mode_hierarchical-SBM_gene-blocks")
# saveRDS(en_head_abs, 'data/enGo_head_abs.Rds')
# saveRDS(en_body_abs, 'data/enGo_body_abs.Rds')

en_head = readRDS('data/enGo_head.Rds')
en_body = readRDS('data/enGo_body.Rds')

en_head_table = ldply(en_head$CP, function(x) x@result) %>% filter(p.adjust < 0.05, Count >= 4)
# write.csv(en_head_table, file  = "en_head.csv")
# en_body_table = ldply(en_body$CP, function(x) x@result) %>% filter(p.adjust < 0.05, Count >= 4)
# write.csv(en_body_table, file  = "en_body.csv")


table_en = table(en_head_abs$summary$n_enrich[en_head_abs$summary$Nested_Level==1]!=0)
table_en
table_en/sum(table_en)

table_en = table(en_head$summary$n_enrich[en_head$summary$Nested_Level==1]!=0)
table_en
table_en/sum(table_en)

table_en = table(en_head$summary$n_enrich[en_head$summary$Nested_Level==2]!=0)
table_en
table_en/sum(table_en)

table_en = table(en_body_abs$summary$n_enrich[en_body_abs$summary$Nested_Level==1]!=0)
table_en

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


wgcna_df =ldply(1:length(en_head$WGCNA), function(id){
      x <- en_head$WGCNA[[id]];
      dplyr::select(x@result, -geneID) %>%
        filter(p.adjust < 0.05, Count >= 4) %>%
        mutate(Module = id-1) %>%
        select(Module, everything())
})
k = table(wgcna_df$ID, wgcna_df$Description)
wgcna_df %>% filter(Module == 2) %>% select(Description, Count)
wgcna_df %>% filter(grepl("synapse", Description)   )
wgcna_df %>% filter(Module==4)

inner_join(en_head_table  %>% filter(grepl("cytoplasmic translation", Description)) %>% select(-geneID) %>% rename(Name=.id),
           en_head$summary, by="Name") %>% filter(Nested_Level == 1) %>% select(Name, ID, Assortatitvity, GeneRatio, p.adjust.x)

ego =
dotplot(ego, showCategory=30)
dotplot(en_head$WGCNA[[5]], showCategory=30)
names = getChild("0-0-0", en_head$summary)[-1]
plot_000 = plot_grid(plotlist = llply(en_head$CP[names], cnetplot, node_label="category",showCategory = 10), ncol = 2, labels = names)
save_plot("g:Dropbox/labbio/articles/NEX_BodyHead_Control-SBM/figures/000_go_map.png", plot_000, base_height = 6, base_asp = 1.1, ncol = 2, nrow = 3)
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
enrichplot::treeplot(enrichplot::pairwise_termsim(en_head$CP$`0-0-0`))

enrichplot::treeplot(enrichplot::pairwise_termsim(en_head$CP$`1-1-1`))
enrichplot::treeplot(enrichplot::pairwise_termsim(en_head$CP$`2-1`))
enrichplot::treeplot(enrichplot::pairwise_termsim(en_head$WGCNA[[2]]))



plot_grid(
emapplot(enrichplot::pairwise_termsim(en_head$WGCNA[[5]])),
empplot(enrichplot::pairwise_termsim(en_head$CP$`0-0-0-0`)))

emplot = function(x) emapplot(enrichplot::pairwise_termsim(x))

install.packages("remotes")
remotes::install_github("GuangchuangYu/enrichplot")

