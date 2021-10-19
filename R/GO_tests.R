source("R/go_functions.R")

# en_head = makeEnrichment("data/output/SBM/clustering/head_weights-spearman_fdr-1e-04_mcmc_mode_hierarchical-SBM_gene-blocks")
# en_body = makeEnrichment("data/output/SBM/clustering/body_weights-spearman_fdr-1e-05_mcmc_mode_hierarchical-SBM_gene-blocks")
# en_head$WGCNA = llply(0:8, getEnrichment,
#                       folder_path="data/output/SBM/clustering/head_weights-spearman_fdr-1e-04_mcmc_mode_hierarchical-SBM_gene-blocks",
#                       clustering = "WGCNA")
# en_body$WGCNA = llply(0:9, getEnrichment,
#                       folder_path="data/output/SBM/clustering/body_weights-spearman_fdr-1e-05_mcmc_mode_hierarchical-SBM_gene-blocks",
#                       clustering = "WGCNA")
# saveRDS(en_head, 'data/enGo_head.Rds')
# saveRDS(en_body, 'data/enGo_body.Rds')

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


goplot_list = llply(en_body$summary$Name[en_body$summary$Nested_Level==3],
                    XGR_plot, en_body$XGR, en_body$summary)
XGR_plot(x="13-1-0", enGo = en_body$XGR, summary = en_body$summary)
#save_plot("go_head_level_11-1-0-super_translation.png", XGR_plot("11-1-0"), base_height = 7, base_asp = 0.25, ncol=4)
goplot_list = llply(en_head$summary$Name[en_head$summary$Nested_Level==3],
                    XGR_plot, en_head$XGR, en_head$summary)

goplot_list = llply(en_body$summary$Name[en_body$summary$Nested_Level==5],
                    XGR_plot, en_body$XGR, en_body$summary)

XGR_plot(x="0-0-0", enGo = en_head$XGR, summary = en_head$summary)



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
plot_000 = plot_grid(plotlist = llply(en_head$CP[names], cnetplot), ncol = 2, labels = names)
save_plot("~/Dropbox/labbio/articles/NEX_BodyHead_Control-SBM/figures/000_go_map.png", plot_000, base_height = 6, base_asp = 1.1, ncol = 2, nrow = 3)
plot_grid(plotlist = llply(en_head$CP[getChild("12-0-0", en_head$summary)[c(-1, -3, -5)]], cnetplot))

getChild("1-1", en_head$summary)[c(-1)]
plot_grid(plotlist = llply(en_head$CP[getChild("1-1", en_head$summary)[c(-1)]], cnetplot))


cnetplot(en_head$CP$`9-4-0`)
cnetplot(en_head$CP$`8-4-0`)

cnetplot(en_head$WGCNA[[7]])
dotplot(en_head$WGCNA[[5]])

?dotplot


rbind(read_csv("data/output/SBM/clustering/head_weights-spearman_fdr-1e-04_mcmc_mode_hierarchical-SBM.csv") %>%
  mutate(tissue = "head"),
read_csv("data/output/SBM/clustering/body_weights-spearman_fdr-1e-05_mcmc_mode_hierarchical-SBM.csv") %>%
  mutate(tissue = "body")) %>% select(Gene, Degree, E_corr, B1:B4, tissue) %>%
  write_csv("block_partition-HeadBody-Control.csv")


bp2 <- clusterProfiler::simplify(en_head$WGCNA[[5]])
cnetplot(clusterProfiler::simplify(en_head$WGCNA[[2]]))
cnetplot(clusterProfiler::simplify(en_head$WGCNA[[3]]))
cnetplot(clusterProfiler::simplify(en_head$WGCNA[[4]]))
cnetplot(clusterProfiler::simplify(en_head$WGCNA[[5]]))
cnetplot(clusterProfiler::simplify(en_head$WGCNA[[6]]))
cnetplot(clusterProfiler::simplify(en_head$WGCNA[[7]]))

