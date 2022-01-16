library(tidyverse)

MMC_body = read_csv("data/output/MMC/mmc.body.output.csv")
MMC_head = read_csv("data/output/MMC/mmc.head.output.csv")

Head_hsbm = read_csv("data/output/SBM/clustering/head_weights-spearman_fdr-1e-04_mcmc_mode_hierarchical-SBM.csv")
Head_hsbm$X1 = NULL
Head_hsbm$tissue = "head"
Head_hsbm = inner_join(Head_hsbm, MMC_head, by= "Gene")


body_hsbm = read_csv("data/output/SBM/clustering/body_weights-spearman_fdr-1e-05_mcmc_mode_hierarchical-SBM.csv")
body_hsbm$X1 = NULL
body_hsbm$tissue = "body"
body_hsbm = inner_join(body_hsbm, MMC_body, by= "Gene")

MMC_HSBM = rbind(body_hsbm, Head_hsbm)

plot = MMC_HSBM %>%
  mutate(B3 = as.factor(B3), B1 = as.factor(B1)) %>%
  ggplot(aes(B2, Module, color = B3)) +
  geom_jitter(width = 0.2, height = 0.2, alpha = 0.3) + facet_wrap(~tissue, scales = "free") +
  scale_y_continuous(breaks = 0:25, labels = c("Not\nclustered", 1:25)) + scale_x_continuous(breaks = 0:25) +
  theme_cowplot() + background_grid() + labs(y = "MMC Modules", x = "SBM\nLevel-2 blocks") +
  scale_color_discrete(name = "SBM\nLevel-3") + theme(legend.position = "bottom")
save_plot("~/Dropbox/labbio/articles/NEX_BodyHead_Control-SBM/figures/MMC_comparison.png", plot, base_height = 5, ncol = 2, base_asp = 1.2)
