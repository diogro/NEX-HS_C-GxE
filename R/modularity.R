source("R/go_functions.R")


en_head = readRDS('data/enGo_head.Rds')
en_body = readRDS('data/enGo_body.Rds')


ggplot(filter(en_head$summary, Nested_Level == 1), aes(Assortatitvity,
                                                       -log2(p.adjust),
                                                       color = Parent)) +
  geom_point() + geom_label_repel(aes(label = Block)) + geom_smooth(color = 1, method = "lm")
ggplot(filter(en_head$summary, Nested_Level == 1), aes(N_genes, -log2(p.adjust), color = Parent)) +
  geom_point() + geom_label_repel(aes(label = Block))
ggplot(filter(en_head$summary, Nested_Level == 1), aes(Assortatitvity, n_enrich, color = Parent)) +
  geom_point() + geom_label_repel(aes(label = Block))
ggplot(filter(en_head$summary, Nested_Level == 1), aes(N_genes, n_enrich, color = Parent)) +
  geom_point() + geom_label_repel(aes(label = Block))
ggplot(filter(en_head$summary, Nested_Level == 1), aes(Assortatitvity, n_enrich_simple, color = Parent)) +
  geom_point() + geom_label_repel(aes(label = Block))
ggplot(filter(en_head$summary, Nested_Level == 1), aes(N_genes, n_enrich_simple, color = Parent)) +
  geom_point() + geom_label_repel(aes(label = Block))




x = en_head$summary %>% filter(Nested_Level == 1) %>% mutate(is_enriched = n_enrich > 0)
x %>%
  mutate(Block = factor(Block, levels = x$Block[order(x$Assortatitvity)])) %>%
ggplot(aes(x=Block, y=Assortatitvity, color = is_enriched)) +
  geom_point(size = 3) + scale_color_colorblind() +
  geom_segment(size = 1, aes(x=Block, xend=Block, y=0, yend=Assortatitvity))

x = en_body$summary %>% filter(Nested_Level == 1) %>% mutate(is_enriched = n_enrich > 0)
x %>%
  mutate(Block = factor(Block, levels = x$Block[order(x$Assortatitvity)])) %>%
  ggplot(aes(x=Block, y=Assortatitvity, color = is_enriched)) +
  geom_point(size = 3) + scale_color_colorblind() +
  geom_segment(size = 1, aes(x=Block, xend=Block, y=0, yend=Assortatitvity))
