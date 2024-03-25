# pak::pkg_install("devtools")
# devtools::install_github("diogro/yamda-r", subdir = "package")
library(yamdar)
if(!require(patchwork)){pak::pkg_install("patchwork"); library(patchwork)}
if(!require(ggplot2)){pak::pkg_install("ggplot2"); library(ggplot2)}
if(!require(cowplot)){pak::pkg_install("cowplot"); library(cowplot)}
if(!require(ggthemes)){pak::pkg_install("ggthemes"); library(ggthemes)}
if(!require(RColorBrewer)){pak::pkg_install("RColorBrewer"); library(RColorBrewer)}
if(!require(evolqg)){pak::pkg_install("evolqg"); library(evolqg)}
if(!require(mvtnorm)){pak::pkg_install("mvtnorm"); library(mvtnorm)}
if(!require(WGCNA)){pak::pkg_install("WGCNA"); library(WGCNA)}

corrPlot = function(M, title = ""){
melted_cormat <- reshape2::melt(M)
    heatmap = ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) +
        geom_tile() +
        scale_fill_gradientn(colours=brewer.pal(11, "RdBu"), 
                             limits = c(-1, 1), 
                              breaks=c(-1, -0.5, 0 , 0.5, 1))  + 
        labs(x = "Traits", y = "") +
        coord_fixed() + 
        theme_tufte() + ggtitle(title) + 
        theme(legend.position = "bottom",
              legend.key.width= unit(1.5, 'cm'), 
              legend.title = element_blank(),
              plot.title = element_text(size = 14),
              axis.title = element_text(size = 12),
              axis.text.y = element_text(size = 6),
              axis.text.x = element_text(size = 6))
    heatmap
}

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

non.mod.cor = mod.cor
between = non.mod.cor[45, 50]
non.mod.cor[1:5, 1:5] = 0.1
diag(non.mod.cor) = 1
write.csv(non.mod.cor, "SBM/non_modular_matrix.csv")


mod.cor[1:10, 1:10]
save_plot("SBM/modular_matrix.pdf", corrPlot(mod.cor))
save_plot("SBM/non_modular_matrix.pdf", corrPlot(non.mod.cor))

## WGCNA

dissTOM = 1 - TOMsimilarity(mod.cor^2)
hierTOM = hclust(as.dist(dissTOM),method="average");
colorDynamicTOM = labels2colors (cutreeDynamic(hierTOM, distM = dissTOM, cutHeight = 0.95, minClusterSize = 3))
pdf("SBM/modular_dendrogram.pdf", width = 3, height = 1.8, pointsize = 6)
plotDendroAndColors(hierTOM,
                    colors=data.frame(colorDynamicTOM),
                    dendroLabels = FALSE, marAll = c(1, 8, 3, 1),
                    main = "")
dev.off()

modular_matrix <- read.csv("SBM/modular_matrix_trimed.csv", header = FALSE)
modular_matrix <- as.matrix(modular_matrix)
dissTOM = 1 - TOMsimilarity(modular_matrix^2)
hierTOM = hclust(as.dist(dissTOM),method="average");
colorDynamicTOM = labels2colors (cutreeDynamic(hierTOM, distM = dissTOM,cutHeight = 0.99, minClusterSize = 3))
pdf("SBM/modular_dendrogram_trimmed.pdf", width = 3, height = 1.8, pointsize = 6)
plotDendroAndColors(hierTOM,
                    colors=data.frame(colorDynamicTOM),
                    dendroLabels = FALSE, marAll = c(1, 8, 3, 1),
                    main = "")
dev.off()

dissTOM = 1 - TOMsimilarity(non.mod.cor^2)
hierTOM = hclust(as.dist(dissTOM),method="average");
colorDynamicTOM = labels2colors (cutreeDynamic(hierTOM, distM = dissTOM,cutHeight = 0.99, minClusterSize = 3))
pdf("SBM/non_modular_dendrogram.pdf", width = 3, height = 1.8, pointsize = 6)
plotDendroAndColors(hierTOM,
                    colors=data.frame(colorDynamicTOM),
                    dendroLabels = FALSE, marAll = c(1, 8, 3, 1))
                    #main = "Gene dendrogram and module colors, TOM dissimilarity")
dev.off()

non_modular_matrix <- read.csv("SBM/non_modular_matrix_trimed.csv", header = FALSE)
non_modular_matrix <- as.matrix(non_modular_matrix)
dissTOM = 1 - TOMsimilarity(non_modular_matrix^2)
hierTOM = hclust(as.dist(dissTOM),method="average");
colorDynamicTOM = labels2colors (cutreeDynamic(hierTOM, distM = dissTOM,cutHeight = 0.99, minClusterSize = 3))
pdf("SBM/non_modular_dendrogram_trimmed.pdf", width = 3, height = 1.8, pointsize = 6)
plotDendroAndColors(hierTOM,
                    colors=data.frame(colorDynamicTOM),
                    dendroLabels = FALSE, marAll = c(1, 8, 3, 1),
                    main = "")
dev.off()
