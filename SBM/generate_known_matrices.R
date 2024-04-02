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


## Figure 1 matrices for cool networks

# matrix with 5 modules

modules = matrix(c(rep(c(1, 0, 0, 0, 0), each = 10),
                   rep(c(0, 1, 0, 0, 0), each = 10),
                   rep(c(0, 0, 1, 0, 0), each = 10),
                   rep(c(0, 0, 0, 1, 0), each = 10),
                   rep(c(0, 0, 0, 0, 1), each = 10)), 50, 5)
modules_z_coef = c(.1, 0.4, .4, .4, .4, 0.4)
mod.cor.5 = calcExpectedMatrix(modules, modules_z_coef)

# Replace values different than 1 with a random value centered on the current value
mod.cor.5[mod.cor.5 != 1] = rnorm(sum(mod.cor.5 != 1), mean = mod.cor.5[mod.cor.5 != 1], sd = 0.1)
# Set a random subsetd of the smaller values with 0
N = sum(mod.cor.5 < 0.3) / 1.1
mod.cor.5[sample(which(mod.cor.5 < 0.3), N)] = 0

N = sum(mod.cor.5 > 0.3 & mod.cor.5 < 0.45) / 2.15
mod.cor.5[sample(which(mod.cor.5 > 0.3 & mod.cor.5 < 0.45), N)] = 0

# Make matrix symmetric
mod.cor.5[lower.tri(mod.cor.5)] = t(mod.cor.5)[lower.tri(mod.cor.5)]
write.csv(mod.cor.5, "SBM/network_architecture/modular_matrix_5modules.csv")

for(m in 1:5){
    index = (m-1)*10 + 1:10
    n_index = index[order(colSums(mod.cor.5[index, index]), decreasing = TRUE)]
    mod.cor.5[index, index] = 
        mod.cor.5[n_index, n_index]
}
save_plot("SBM/network_architecture/modular_matrix_5modules.png", corrPlot(mod.cor.5) + theme_void() + theme(legend.position = "none"))

# Core-periphery matrix

modules = matrix(c(rep(c(1, 0, 0, 0, 0), each = 10),
                   rep(c(1, 1, 1, 1, 1), each = 10),
                   rep(c(0, 1, 1, 1, 1), each = 10),
                   rep(c(0, 0, 1, 1, 1), each = 10),
                   rep(c(0, 0, 0, 1, 1), each = 10),
                   rep(c(0, 0, 0, 0, 1), each = 10)), 50, 6)
modules_z_coef = c(.45, 
                   .2, -.2, -.1, -0.1, -0.05, -0.05)
mod.cor.cp = calcExpectedMatrix(modules, modules_z_coef)

# Replace values different than 1 with a random value centered on the current value
mod.cor.cp[mod.cor.cp != 1] = rnorm(sum(mod.cor.cp != 1), mean = mod.cor.cp[mod.cor.cp != 1], sd = 0.1)
# Set a random subsetd of the smaller values with 0
N = sum(mod.cor.cp < 0.3) / 1.2
mod.cor.cp[sample(which(mod.cor.cp < 0.3), N)] = 0

# N = sum(mod.cor.cp > 0.35) / 2.5
# mod.cor.cp[sample(which(mod.cor.cp > 0.35), N)] = 0

# Make matrix symmetric
mod.cor.cp[lower.tri(mod.cor.cp)] = t(mod.cor.cp)[lower.tri(mod.cor.cp)]
write.csv(mod.cor.cp, "SBM/network_architecture/cp_matrix.csv")

for(m in 1:5){
    index = (m-1)*10 + 1:10
    n_index = index[order(colSums(mod.cor.cp[index, index]), decreasing = TRUE)]
    mod.cor.cp[index, index] = 
        mod.cor.cp[n_index, n_index]
}
diag(mod.cor.cp) = 1
save_plot("SBM/network_architecture/cp_matrix.png", corrPlot(mod.cor.cp) + theme_void() + theme(legend.position = "none"))



# Disassortative matrix

modules = matrix(c(rep(c(1, 0, 0, 0, 0), each = 10),
                   rep(c(0, 1, 0, 0, 0), each = 10),
                   rep(c(0, 0, 1, 0, 0), each = 10),
                   rep(c(0, 0, 0, 1, 0), each = 10),
                   rep(c(0, 0, 0, 0, 1), each = 10)), 50, 6)
modules_z_coef = c(.5, 
                   -.5, -.5, -.5, -0.5, -0.5)
mod.cor.dis = calcExpectedMatrix(modules, modules_z_coef)

for(i in 1:20){
    for(j in 21:50){
        mod.cor.dis[i, j] = 0
    }
}

for(i in 31:40){
    for(j in 41:50){
        mod.cor.dis[i, j] = 0
    }
}

for(i in 11:20){
    for(j in 41:50){
        mod.cor.dis[i, j] = 0.46
    }
}

# Replace values different than 1 with a random value centered on the current value
mod.cor.dis[mod.cor.dis != 1] = rnorm(sum(mod.cor.dis != 1), mean = mod.cor.dis[mod.cor.dis != 1], sd = 0.1)
# Set a random subsetd of the smaller values with 0
N = sum(mod.cor.dis < 0.3) / 1.2
mod.cor.dis[sample(which(mod.cor.dis < 0.3), N)] = 0

 N = sum(mod.cor.dis > 0.35) / 2.5
 mod.cor.dis[sample(which(mod.cor.dis > 0.35), N)] = 0

# Make matrix symmetric
mod.cor.dis[lower.tri(mod.cor.dis)] = t(mod.cor.dis)[lower.tri(mod.cor.dis)]
write.csv(mod.cor.dis, "SBM/network_architecture/dis_matrix.csv")

for(m in 1:5){
    index = (m-1)*10 + 1:10
    n_index = index[order(colSums(mod.cor.dis[index, index]), decreasing = TRUE)]
    mod.cor.dis[index, index] = 
        mod.cor.dis[n_index, n_index]
}
diag(mod.cor.dis) = 1
save_plot("SBM/network_architecture/dis_matrix.png", corrPlot(mod.cor.dis) + theme_void() + theme(legend.position = "none"))


























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


body_matrix <- read.csv("SBM/body_matrix_trimed.csv", header = FALSE)
body_matrix <- as.matrix(body_matrix)
dissTOM = 1 - TOMsimilarity(body_matrix^2)
hierTOM = hclust(as.dist(dissTOM),method="average");
colorDynamicTOM = labels2colors (cutreeDynamic(hierTOM, distM = dissTOM,cutHeight = 0.99, minClusterSize = 3))
pdf("SBM/body_matrix_trimed_WGCNA.pdf", width = 3, height = 1.8, pointsize = 6)
plotDendroAndColors(hierTOM,
                    colors=data.frame(colorDynamicTOM),
                    dendroLabels = FALSE, marAll = c(1, 8, 3, 1),
                    main = "")
dev.off()

head_matrix <- read.csv("SBM/head_matrix_trimed.csv", header = FALSE)
head_matrix <- as.matrix(head_matrix)
dissTOM = 1 - TOMsimilarity(head_matrix^2)
hierTOM = hclust(as.dist(dissTOM),method="average");
colorDynamicTOM = labels2colors (cutreeDynamic(hierTOM, distM = dissTOM,cutHeight = 0.99, minClusterSize = 3))
pdf("SBM/head_matrix_trimed_WGCNA.pdf", width = 3, height = 1.8, pointsize = 6)
plotDendroAndColors(hierTOM,
                    colors=data.frame(colorDynamicTOM),
                    dendroLabels = FALSE, marAll = c(1, 8, 3, 1),
                    main = "")
dev.off()