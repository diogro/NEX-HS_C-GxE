source(here::here("R/go_functions.R"))

#plots_path = "B:/Dropbox/labbio/articles/NEX_BodyHead_Control-SBM/figures/"
plots_path = "~/Dropbox/labbio/articles/NEX_BodyHead_Control-SBM/figures/"


# en_head = makeEnrichment("data/output/SBM/clustering/head_weights-spearman_fdr-1e-02_mcmc_mode_hierarchical-SBM_gene-blocks")
# en_body = makeEnrichment("data/output/SBM/clustering/body_weights-spearman_fdr-1e-03_mcmc_mode_hierarchical-SBM_gene-blocks")
# en_head$WGCNA = llply(0:7, getEnrichment,
#                       folder_path="data/output/SBM/clustering/head_weights-spearman_fdr-1e-02_mcmc_mode_hierarchical-SBM_gene-blocks",
#                        clustering = "WGCNA")
# en_body$WGCNA = llply(0:7, getEnrichment,
#                       folder_path="data/output/SBM/clustering/body_weights-spearman_fdr-1e-03_mcmc_mode_hierarchical-SBM_gene-blocks",
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

en_body$summary %>%
  filter(Nested_Level == 1) %>%
  arrange(desc(Assortatitvity)) %>% head()

CP_print("4-0-0-0", en_head$CP, en_head$summary) %>%
  filter(Block == "4-0-0-0") 
CP_print("8-7-2-0-0", en_body$CP, en_body$summary)

#%%
### Table 1: Go enrichment of body and head

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
#%%

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

ggplot(en_body$summary, aes(N_genes, n_enrich)) + geom_point() + theme_cowplot() + facet_wrap(~Nested_Level, scales="free")

wgcna_df_head %>%
  group_by(Module) %>%
  filter(Module > 0) %>%
  count()
wgcna_df_body %>%
  group_by(Module) %>%
  filter(Module > 0) %>%
  count()

wgcna_df_body %>% filter(grepl("cytoplasmic translation", Description))
wgcna_df_head %>% filter(grepl("ATP", Description))


WGCNA_head_all = llply(en_head$WGCNA[-1], barplot)
plot_grid(plotlist = WGCNA_head_all, ncol = 3, labels = 1:7, scale = 0.9)

WGCNA_body_all = llply(en_body$WGCNA[-1], barplot)
plot_grid(plotlist = WGCNA_body_all, ncol = 3, labels = 1:7, scale = 0.9)

names = en_head$summary %>% filter(Nested_Level == 3)
Level_3_list = llply(en_head$CP[names$Name], barplot)
plot_grid(plotlist = Level_3_list, ncol = 3, labels = names$Name, scale = 0.9)

names = en_head$summary %>% filter(Nested_Level == 2, n_enrich > 0)
Level_2_list = llply(en_head$CP[names$Name], barplot)
plot_grid(plotlist = Level_2_list, ncol = 4, labels = names$Name, scale = 0.9)

names = en_body$summary %>% filter(Nested_Level == 3)
Level_3_list = llply(en_body$CP[names$Name], barplot)
plot_grid(plotlist = Level_3_list, ncol = 3, labels = names$Name, scale = 0.9)

names = en_body$summary %>% filter(Nested_Level == 2, n_enrich > 0)
Level_2_list = llply(en_body$CP[names$Name], barplot)
plot_grid(plotlist = Level_2_list, ncol = 4, labels = names$Name, scale = 0.9)

wgcna_df_head %>% filter(Module==2)

wgcna_df_body %>%
  group_by(Module) %>%
  count()
wgcna_df_head %>%
  group_by(Module) %>%
  count()

names = c(getChild("6-0-0-0", en_head$summary |> filter(n_enrich > 0))[-1],
          getChild("4-0-0-0", en_head$summary |> filter(n_enrich > 0))[-1],
          getChild("0-0-0-0", en_head$summary |> filter(n_enrich > 0))[-1])[c(-4, -6)]
SBM_neuro_list = llply(en_head$CP[names], barplot)
WGCNA_neuro = barplot(en_head$WGCNA[[5]])

neuro_list = SBM_neuro_list; neuro_list[[length(neuro_list)+1]] = WGCNA_neuro
plot_000 = plot_grid(plotlist = neuro_list, ncol = 2, labels = c(names, "WGCNA module 4"), scale = 0.95)
plot_000
save_plot(file.path(plots_path, "000_go_map.png"), plot_000, base_height = 4.3,
          base_asp = 1.8, ncol = 2, nrow = 4)

names = getChild("8-4-1-1", en_head$summary |> filter(n_enrich > 0))[-1]
SBM_translation_list = llply(en_head$CP[names], barplot)
WGCNA_translation = barplot(en_head$WGCNA[[3]])

translation_list = SBM_translation_list; translation_list[[length(translation_list)+1]] = WGCNA_translation
plot_411 = plot_grid(plotlist = translation_list, ncol = 2, labels = c(names, "WGCNA module 4"), scale = 0.9)
save_plot(file.path(plots_path, "411_translation_map.png"), plot_411, base_height = 5,
          base_asp = 1.5, ncol = 2, nrow = 3)

translate_head = inner_join(en_head_table  %>% filter(grepl("cytoplasmic translation", Description) |
                                                        grepl("^translation$", Description)) %>%
                              select(-geneID, -pvalue),
           en_head$summary, by="Name") %>% filter(Nested_Level == 1)
translate_body = inner_join(en_body_table  %>% filter(grepl("cytoplasmic translation", Description) |
                                                        grepl("^translation$", Description)) %>%
                              select(-geneID, -pvalue),
                            en_body$summary, by="Name") %>% filter(Nested_Level == 1)
translate_all = rbind(translate_body, translate_head) %>%
  filter(N_genes < 100) 
summary(translate_all$Assortatitvity)

neuro_head = inner_join(en_head_table  %>% filter(grepl("G protein-coupled receptor", Description)) %>%
                              select(-geneID, -pvalue),
                            en_head$summary, by="Name") %>% filter(Nested_Level == 1)

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

mean(c(translate_head$Assortativity, translate_body$Assortativity))
min( c(translate_head$Assortativity, translate_body$Assortativity))
max( c(translate_head$Assortativity, translate_body$Assortativity))


names_head = translate_head$Name
names_body = translate_body$Name
SBM_translation_list = llply(en_head$CP[names_head], cnetplot, node_label="category",showCategory = 10)
SBM_translation_list = append(SBM_translation_list, llply(en_head$CP[names_body], cnetplot, node_label="category",showCategory = 10))
plot_translation= plot_grid(plotlist = SBM_translation_list, ncol = 2,
                            labels = c(paste("Head -", names_head), paste("Body -", names_head)), scale = 0.9)
save_plot(file.path(plots_path, "Translation_go_map.png"), plot_translation, base_height = 5,
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

# Some non-assortative modules in the head:

head_non_assort = en_head$summary %>% 
  filter(Assortativity < 0 & n_enrich > 0)

en_head_table %>%
  filter(Name %in% head_non_assort$Name)

pl = llply(head_non_assort$Name, auto_barplot, en_head$CP)
p = plot_grid(plotlist = pl, labels = head_non_assort$Name)
save_plot("head_non_assortative_options.png", p, ncol = 2, nrow = 2, base_height = 5)

body_non_assort = en_body$summary %>% 
  filter(Assortativity < -0.006 & n_enrich > 0 & Nested_Level == 1 & N_genes < 30)

pl = llply(body_non_assort$Name, auto_barplot, en_body$CP)
p = plot_grid(plotlist = pl, labels = body_non_assort$Name)
save_plot("body_non_assortative_options.png", p, ncol = 3, nrow = 2, base_height = 5)



### Hypergeometric test Assortativity and Cytoplasmic Translation
# install.packages("ggpubr")
library(ggpubr)

# Body
{
translation_assortativity = inner_join(en_body_table  %>% 
           mutate(Translation = grepl("cytoplasmic translation", Description) | grepl("^translation$", Description)) %>%
           select(-geneID, -pvalue),
           en_body$summary, by="Name") %>% 
           filter(Nested_Level == 1) %>% 
           select(Name, Assortativity, Translation) %>% 
           unique()

compare_means(Assortativity ~ Translation, data = translation_assortativity, method = "wilcox.test")

my_comparisons <- list(c(1, 2))
p <- ggboxplot(translation_assortativity, x = "Translation", y = "Assortativity",
                    add = "jitter") + scale_x_discrete(labels = c("Other terms", "Cytoplasmic translation")) +
                    labs(x = "GO annotation")
#  Add p-value
p_body = p + 
  stat_compare_means(label.y = 0.036, label.x = 0.7) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test") +
  ggtitle("A. Body")


# Head
translation_assortativity = inner_join(en_head_table  %>% 
           mutate(Translation = grepl("cytoplasmic translation", Description) | grepl("^translation$", Description)) %>%
           select(-geneID, -pvalue),
           en_head$summary, by="Name") %>% 
           filter(Nested_Level == 1) %>% 
           select(Name, Assortativity, Translation) %>% 
           unique()

compare_means(Assortativity ~ Translation, data = translation_assortativity, method = "wilcox.test")

my_comparisons <- list(c(1, 2))
p <- ggboxplot(translation_assortativity, x = "Translation", y = "Assortativity",
                    add = "jitter") + scale_x_discrete(labels = c("Other terms", "Cytoplasmic translation")) +
                    labs(x = "GO annotation")
#  Add p-value
p_head = p + 
  stat_compare_means(label.y = 0.17, label.x = 0.7) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test") +
  ggtitle("B. Head")
# Change method
#p + stat_compare_means(method = "t.test")
panel = p_body + p_head 
save_plot("assortativity_cytoplasmic_translation.png", panel, base_height = 4, ncol = 2, base_asp = 1)
}

{
translation_assortativity = inner_join(en_body_table  %>% 
           mutate(Respiration = grepl("respiration", Description) | 
                    grepl("ATP", Description) | 
                    grepl("oxidative", Description) | 
                    grepl("carboxylic", Description)) %>%
           select(-geneID, -pvalue),
           en_body$summary, by="Name") %>% 
           filter(Nested_Level == 1) %>% 
           select(Name, Assortativity, Respiration) %>% 
           unique()

compare_means(Assortativity ~ Respiration, data = translation_assortativity, method = "wilcox.test")

my_comparisons <- list(c(1, 2))
p <- ggboxplot(translation_assortativity, x = "Respiration", y = "Assortativity",
                    add = "jitter") + scale_x_discrete(labels = c("Other terms", "Cell Respiration")) +
                    labs(x = "GO annotation")
#  Add p-value
p_body = p + 
  stat_compare_means(label.y = 0.036, label.x = 0.7) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test")


inner_join(en_body_table  %>% 
           filter((grepl("respiration", Description) | 
                    grepl("ATP", Description) | 
                    grepl("oxidative", Description) | 
                    grepl("carboxylic", Description))) %>%
           select(-geneID, -pvalue),
           en_body$summary, by="Name") %>% 
           filter(Nested_Level == 1) %>% 
           select(Name, Description) 

# Head
translation_assortativity = inner_join(en_head_table  %>% 
           mutate(Respiration = grepl("respiration", Description) | 
                    grepl("ATP", Description) | 
                    grepl("oxidative", Description) | 
                    grepl("carboxylic", Description)) %>%
           select(-geneID, -pvalue),
           en_head$summary, by="Name") %>% 
           filter(Nested_Level == 1) %>% 
           select(Name, Assortativity, Respiration) %>% 
           unique()

compare_means(Assortativity ~ Respiration, data = translation_assortativity, method = "wilcox.test")

my_comparisons <- list(c(1, 2))
p <- ggboxplot(translation_assortativity, x = "Respiration", y = "Assortativity",
                    add = "jitter") + scale_x_discrete(labels = c("Other terms", "Cell Respiration")) +
                    labs(x = "GO annotation")
#  Add p-value
p_head = p + 
  stat_compare_means(label.y = 0.17, label.x = 0.7) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test")
# Change method
panel = p_body + p_head 
save_plot("assortativity_Respiration.png", panel, base_height = 4, ncol = 2, base_asp = 1)
}
