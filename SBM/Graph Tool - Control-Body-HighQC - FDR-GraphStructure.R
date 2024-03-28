library(tidyverse)
library(cowplot)
library(patchwork)

data <- read_csv("data/output/SBM/FDRCutOff/body_cutOff_Ngenes_density.csv")
p1 = data %>%
  ggplot(aes(x = (FDRcutOff), y = N)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 1e-3) +
  scale_y_continuous(limits = c(0, 6000)) + 
  labs(title = "Body - FDR vs Ngenes",
       x = "FDR",
       y = "Ngenes",
       color = "Graph Structure") +
  theme(legend.position = "bottom") +  
  theme_cowplot()

p2 = data %>%
  ggplot(aes(x = (FDRcutOff), y = density)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 1e-3) + 
  scale_y_continuous(limits = c(0, 1)) + 
  labs(title = "Body - FDR vs Density",
       x = "FDR",
       y = "Density",
       color = "Graph Structure") +
  theme(legend.position = "bottom") +  
  theme_cowplot()

data <- read_csv("data/output/SBM/FDRCutOff/head_cutOff_Ngenes_density.csv")

p3 = data %>%
  ggplot(aes(x = (FDRcutOff), y = N)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 1e-2) + 
  scale_y_continuous(limits = c(0, 6000)) + 
    labs(title = "Head - FDR vs Ngenes",
       x = "FDR",
       y = "Ngenes",
       color = "Graph Structure") +
  theme(legend.position = "bottom") +  
  theme_cowplot()

p4 = data %>%
  ggplot(aes(x = (FDRcutOff), y = density)) +
    geom_point() +
    geom_line() +
    geom_vline(xintercept = 1e-2) + 
    scale_y_continuous(limits = c(0, 1)) + 
    labs(title = "Head - FDR vs Density",
         x = "FDR",
         y = "Density",
         color = "Graph Structure") +
    theme(legend.position = "bottom") +  
    theme_cowplot()

panel = plot_grid(p1, p3, p2, p4, ncol = 2)
save_plot("data/output/SBM/FDRCutOff/cutOff_Ngenes.png", panel, base_height = 10, base_width = 10)

