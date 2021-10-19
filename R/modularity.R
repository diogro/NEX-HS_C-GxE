source("R/go_functions.R")


en_head = readRDS('data/enGo_head.Rds')
en_body = readRDS('data/enGo_body.Rds')


ggplot(filter(en_head$summary, Nested_Level == 1), aes(Assortatitvity,
                                                       -log2(p.adjust),
                                                       color = Parent)) +
  geom_point() + geom_label_repel(aes(label = Block))
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



# Create data
data <- data.frame(x=seq(1,30), y=abs(rnorm(30)))

# Plot

x = en_head$summary %>%
  filter(Nested_Level == 1)
x %>%
  mutate(Block = factor(Block, levels = x$Block[order(x$Assortatitvity)])) %>%
ggplot(aes(x=Block, y=Assortatitvity, color = -log10(p.adjust))) +
  geom_point() +
  geom_segment( aes(x=Block, xend=Block, y=0, yend=Assortatitvity))

