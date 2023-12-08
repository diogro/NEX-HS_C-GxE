source("R/go_functions.R")


en_head = readRDS('data/enGo_head_fdr-1e-02.Rds')
en_body = readRDS('data/enGo_body_fdr-1e-03.Rds')


ggplot(filter(en_head$summary, Nested_Level == 1), aes(Assortativity,
                                                       -log2(p.adjust),
                                                       color = Parent)) +
  geom_point() + geom_label_repel(aes(label = Block)) + geom_smooth(color = 1, method = "lm")
ggplot(filter(en_head$summary, Nested_Level == 1), aes(N_genes, -log2(p.adjust), color = Parent)) +
  geom_point() + geom_label_repel(aes(label = Block))
ggplot(filter(en_head$summary, Nested_Level == 1), aes(Assortativity, n_enrich, color = Parent)) +
  geom_point() + geom_label_repel(aes(label = Block))
ggplot(filter(en_head$summary, Nested_Level == 1), aes(N_genes, n_enrich, color = Parent)) +
  geom_point() + geom_label_repel(aes(label = Block))
ggplot(filter(en_head$summary, Nested_Level == 1), aes(Assortativity, n_enrich_simple, color = Parent)) +
  geom_point() + geom_label_repel(aes(label = Block))
ggplot(filter(en_head$summary, Nested_Level == 1), aes(N_genes, n_enrich_simple, color = Parent)) +
  geom_point() + geom_label_repel(aes(label = Block))


pak::pkg_install("ggthemes")
library(ggthemes)
{
  x = en_head$summary %>% filter(Nested_Level == 1) %>% mutate(is_enriched = n_enrich > 0)
  p_head = x %>%
    mutate(Block = factor(Block, levels = x$Block[order(x$Assortativity)])) %>%
    ggplot(aes(x=Block, y=Assortativity, color = is_enriched)) +
      geom_point(size = 1) + scale_color_colorblind(name = "GO Enriched") +
      scale_x_discrete(guide = guide_axis(n.dodge = 2)) + 
      geom_segment(linewidth = 0.5, aes(x=Block, xend=Block, y=0, yend=Assortativity)) +
      theme_cowplot() + geom_hline(yintercept = 0) + theme(legend.position = "none") +
      theme(axis.text.x = element_text(angle=45, vjust=0.5, size = 4.5, margin = margin(t = 0, r = 0, b = 0.2, l = 0) ),
            plot.title = element_text(size = 12), 
            axis.text.y = element_text(size = 8),
            axis.title = element_text(size = 8),
            legend.text = element_text(size = 6),
            legend.title = element_text(size = 6),
            axis.ticks.x = element_line(size = .3),
            axis.ticks.length=unit(.07, "cm")) +
      labs(x = "Level-1 Blocks", y = "Assortativity")


  x = en_body$summary %>% filter(Nested_Level == 1) %>% mutate(is_enriched = n_enrich > 0)
  p_body = x %>%
    mutate(Block = factor(Block, levels = x$Block[order(x$Assortativity)])) %>%
    ggplot(aes(x=Block, y=Assortativity, color = is_enriched)) +
      geom_point(size = 1) + scale_color_colorblind(name = "GO Enriched", labels = c("No", "Yes")) +
      geom_segment(linewidth = 0.5, aes(x=Block, xend=Block, y=0, yend=Assortativity)) +
      scale_x_discrete(guide = guide_axis(n.dodge = 2)) + 
      theme_cowplot() + geom_hline(yintercept = 0) +
      theme(legend.position = c(0.7, 0.9)) +
      theme(axis.text.x = element_text(angle=45, vjust=0.5, size = 4.5, margin = margin(t = 0, r = 0, b = 0.2, l = 0) ),
            plot.title = element_text(size = 12), 
            axis.text.y = element_text(size = 8),
            axis.title = element_text(size = 8),
            legend.text = element_text(size = 6),
            legend.title = element_text(size = 6),
            axis.ticks.x = element_line(size = .3),
            axis.ticks.length=unit(.07, "cm")) +
      guides(colour = guide_legend(override.aes = list(size=2)))
      labs(x = "Level-1 Blocks", y = "Assortativity")


  x = en_head$summary %>%
    group_by(Nested_Level) %>%
    summarise(mean(Assortativity))
  names(x) = c("Level", "Modularity")
  mod_head = ggplot(x, aes(Level, Modularity)) + geom_point(size = 0.5) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) + 
    geom_line(linewidth = 0.5) + theme_cowplot() + background_grid() + labs(x = "Nested Block Level") + 
    theme(axis.text = element_text(size = 5, margin = margin(t = 0, r = 0, b = 0., l = 0)), 
          axis.title = element_text(size = 5, margin = margin(t = 0, r = 0, b = 0., l = 0)), 
          axis.ticks.x = element_line(size = .3),
          axis.ticks.length=unit(.05, "cm"))

  x = en_body$summary %>%
    group_by(Nested_Level) %>%
    summarise(mean(Assortativity))
  names(x) = c("Level", "Modularity")
  mod_body = ggplot(x, aes(Level, Modularity)) + geom_point(size = 0.5) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) + 
    geom_line(linewidth = 0.5) + theme_cowplot() + background_grid() + labs(x = "Nested Block Level") + 
    theme(axis.text = element_text(size = 5), axis.title = element_text(size = 5), axis.ticks = element_blank())


  pb = ggdraw(p_body) +
            draw_plot(mod_body, .115, .5, .5, .5)
  ph = ggdraw(p_head) +
    draw_plot(mod_head, .115, .4, .5, .5)

  plot = plot_grid(pb, ph, ncol  = 1, labels = c("A. Body", "B. Head"), scale = 0.9, label_size = 8)
  plot
  save_plot("~/Dropbox/labbio/articles/SBM_manuscript/figures//assortativity.png",
            plot, base_width = 5.2, base_height = 3.5)
}
