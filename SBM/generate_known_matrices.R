library(cowplot)
library(yamdar)
if(!require(patchwork)){pak::pkg_install("patchwork"); library(patchwork)}
if(!require(ggplot2)){pak::pkg_install("ggplot2"); library(ggplot2)}
if(!require(cowplot)){pak::pkg_install("cowplot"); library(cowplot)}
if(!require(ggthemes)){pak::pkg_install("ggthemes"); library(ggthemes)}
if(!require(RColorBrewer)){pak::pkg_install("RColorBrewer"); library(RColorBrewer)}
library(evolqg)
library(mvtnorm)
pak::pkg_install("WGCNA")
library(WGCNA)

modules = matrix(c(rep(c(1, 0, 0, 0, 0), each = 10),
                   rep(c(0, 1, 0, 0, 0), each = 10),
                   rep(c(0, 0, 1, 0, 0), each = 10),
                   rep(c(0, 0, 0, 1, 0), each = 10),
                   rep(c(0, 0, 0, 0, 1), each = 10),
                   c(rep(1, 20), rep(0, 30)),
                   c(rep(0, 20), rep(1, 20), rep(0, 10))), 50, 7)
modules_z_coef = c(.1, 0.2, .2, .2, .2, .2, .3, 0.3)
mod.cor = calcExpectedMatrix(modules, modules_z_coef)
write.csv(mod.cor, "SBM/modular_matrix.csv")
png("test.png")
corrPlot(mod.cor)
dev.off()

between = mod.cor[45, 50]
mod.cor[1:5, 1:5] = 0.1
diag(mod.cor) = 1
write.csv(mod.cor, "SBM/non_modular_matrix.csv")
pop = rmvnorm(1000, sigma = mod.cor)


mod.cor[1:10, 1:10]
png("test.png")
corrPlot(mod.cor)
dev.off()
corrPlot = function(M, title = ""){
melted_cormat <- reshape2::melt(M)
    heatmap = ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) +
        geom_tile() +
        scale_fill_gradientn(colours=brewer.pal(11, "RdBu"), 
                             limits = c(-1, 1), 
                              breaks=c(-1, -0.5, 0 , 0.5, 1))  + 
        labs(x = "Traits", y = "") +
        theme_tufte() + ggtitle(title) + 
        theme(legend.position = "bottom",
              legend.key.width= unit(1.5, 'cm'), 
              legend.title = element_blank(),
              plot.title = element_text(size = 14),
              axis.title = element_text(size = 12),
              axis.text.y = element_text(size = 6),
              axis.text.x = element_text(size = 6, angle = 90))
    heatmap
}

hclust(pop)
