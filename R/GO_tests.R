source("R/go_functions.R")

en_head = makeEnrichment("data/output/SBM/clustering/head_weights-spearman_fdr-1e-04_mcmc_mode_hierarchical-SBM_gene-blocks")
en_body = makeEnrichment("data/output/SBM/clustering/body_weights-spearman_fdr-1e-05_mcmc_mode_hierarchical-SBM_gene-blocks")
saveRDS(en_head, 'data/enGo_head.Rds')
saveRDS(en_body, 'data/enGo_body.Rds')

en_head_abs = makeEnrichment("data/output/SBM/clustering/head_weights-spearman_fdr-1e-04_absolute_mcmc_mode_hierarchical-SBM_gene-blocks")
en_body_abs = makeEnrichment("data/output/SBM/clustering/body_weights-spearman_fdr-1e-05_absolute_mcmc_mode_hierarchical-SBM_gene-blocks")
saveRDS(en_head_abs, 'data/enGo_head_abs.Rds')
saveRDS(en_body_abs, 'data/enGo_body_abs.Rds')

en_head = readRDS('data/enGo_head.Rds')
en_body = readRDS('data/enGo_body.Rds')

en_head_table = ldply(en_head$CP, function(x) x@result) %>% filter(p.adjust < 0.05, Count >= 4)
write.csv(en_head_table, file  = "en_head.csv")
en_body_table = ldply(en_body$CP, function(x) x@result) %>% filter(p.adjust < 0.05, Count >= 4)
write.csv(en_body_table, file  = "en_body.csv")



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

ggplot(filter(block_summary, Nested_Level == 1), aes(Assortatitvity,
                                                     -log2(p.adjust),
                                                     color = Parent)) +
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

CP_print("9-4-0", enGo = en_head$CP, en_head$summary)
en_head$XGR$`1-1-1`
XGR_plot(x = "12-0-0", enGo = en_head$XGR_CC, summary = en_head$summary)
XGR_plot(x = "0-0-0", enGo = en_head$XGR_CC, summary = en_head$summary)
save_plot("go_head_level_0-0-0-vision-neuro-signaling.png", XGR_plot("0-0-0", en_head$XGR, en_head$summary),
          base_height = 7, base_asp = 0.25, ncol=5)
XGR_plot(x = "0-0-0", enGo = en_head$XGR_CC, summary = en_head$summary)
save_plot("go_head_level_0-0-high-level-neuro.png", XGR_plot("0-0", en_head$XGR, en_head$summary),
          base_height = 7, base_asp = 0.25, ncol=5)
save_plot("go_head_level_12-0-0-nueromuscular-signaling.png", XGR_plot("12-0-0", en_head$XGR, en_head$summary),
          base_height = 7, base_asp = 0.25, ncol=5)
save_plot("go_head_level_12-0-0-nueromuscular-signaling_cc.png", XGR_plot("12-0-0", en_head$XGR_CC, en_head$summary, "CC"),
          base_height = 7, base_asp = 0.25, ncol=5)

XGR_plot("0-0-0", en_head$XGR_CC, en_head$summary)
save_plot("go_head_level_0-0-0-vision-neuro-signaling-CC.png", XGR_plot("0-0-0", en_head$XGR_CC, en_head$summary, "CC"), base_height = 7, base_asp = 0.25, ncol=5)


XGR_plot(x = "2-2-1", enGo = en_head$XGR, summary = en_head$summary)
XGR_plot(x = "1-1-1", enGo = en_head$XGR_CC, summary = en_head$summary)
save_plot("go_head_level_1-1-1-super_translation.png", XGR_plot("1-1-1", en_head$XGR, en_head$summary),
          base_height = 7, base_asp = 0.25, ncol=5)
save_plot("go_head_level_1-1-1-super_translation_CC.png", XGR_plot("1-1-1", en_head$XGR_CC, en_head$summary, "CC"),
          base_height = 7, base_asp = 0.25, ncol=5)

save_plot("go_head_level_6-1-1-less_translation.png", XGR_plot("6-1-1", en_head$XGR, en_head$summary),
          base_height = 7, base_asp = 0.25, ncol=5)
XGR_plot(x = "6-1-1", enGo = en_head$XGR_CC, summary = en_head$summary, "CC")

XGR_plot(x = "4-0", enGo = en_head$XGR, summary = en_head$summary)
XGR_plot(x = "4-0", enGo = en_head$XGR_CC, summary = en_head$summary)
save_plot("go_head_level_4-0-muscle-Krebs.png", XGR_plot("4-0", en_head$XGR, en_head$summary),
          base_height = 7, base_asp = 0.25, ncol=5)
save_plot("go_head_level_9-4-0-muscle-Krebs.png", XGR_plot("9-4-0", en_head$XGR, en_head$summary, fdr = 0.1),
          base_height = 7, base_asp = 0.25, ncol=5)
save_plot("go_head_level_1-1-1-super_translation_CC.png", XGR_plot("1-1-1", en_head$XGR_CC, en_head$summary, "CC"),
          base_height = 7, base_asp = 0.25, ncol=5)
