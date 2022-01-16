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


library(ggthemes)

x = en_head$summary %>% filter(Nested_Level == 1) %>% mutate(is_enriched = n_enrich > 0)
p_head = x %>%
  mutate(Block = factor(Block, levels = x$Block[order(x$Assortatitvity)])) %>%
  ggplot(aes(x=Block, y=Assortatitvity, color = is_enriched)) +
    geom_point(size = 3) + scale_color_colorblind(name = "GO Enriched") +
    geom_segment(size = 1, aes(x=Block, xend=Block, y=0, yend=Assortatitvity)) +
    theme_cowplot() + geom_hline(yintercept = 0) + theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle=45, vjust=0.6)) +
    labs(x = "Level 1 Block", y = "Assortativity")


x = en_body$summary %>% filter(Nested_Level == 1) %>% mutate(is_enriched = n_enrich > 0)
p_body = x %>%
  mutate(Block = factor(Block, levels = x$Block[order(x$Assortatitvity)])) %>%
  ggplot(aes(x=Block, y=Assortatitvity, color = is_enriched)) +
    geom_point(size = 3) + scale_color_colorblind(name = "GO Enriched") +
    geom_segment(size = 1, aes(x=Block, xend=Block, y=0, yend=Assortatitvity)) +
    theme_cowplot() + geom_hline(yintercept = 0) +
    theme(legend.position = c(0.7, 0.9)) +
    theme(axis.text.x = element_text(angle=45, vjust=0.6)) +
    labs(x = "Level 1 Block", y = "Assortativity")


x = en_head$summary %>%
  group_by(Nested_Level) %>%
  summarise(mean(Assortatitvity))
names(x) = c("Level", "Modularity")
mod_head = ggplot(x, aes(Level, Modularity)) + geom_point() +
  geom_line() + theme_cowplot() + background_grid() + labs(x = "Nested Block Level")

x = en_body$summary %>%
  group_by(Nested_Level) %>%
  summarise(mean(Assortatitvity))
names(x) = c("Level", "Modularity")
mod_body = ggplot(x, aes(Level, Modularity)) + geom_point() +
  geom_line() + theme_cowplot() + background_grid() + labs(x = "Nested Block Level")


pb = ggdraw(p_body) +
            draw_plot(mod_body, .1, .5, .4, .4)
ph = ggdraw(p_head) +
  draw_plot(mod_head, .1, .4, .4, .4)

plot = plot_grid(pb, ph, ncol  = 1, labels = c("A. Body", "B. Head"), scale = 0.9)
plot
save_plot("~/Dropbox/labbio/articles/NEX_BodyHead_Control-SBM/figures//assortativity.png",
          plot, base_height = 5, ncol = 1, nrow = 2, base_asp = 3.5, )
