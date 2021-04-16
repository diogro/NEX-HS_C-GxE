if(!require(WGCNA)){BiocManager::install("WGCNA"); library(WGCNA)}
if(!require(plyr)){install.packages("plyr"); library(plyr)}
if(!require(tidyverse)){install.packages("tidyverse"); library(tidyverse)}
if(!require(superheat)){install.packages("superheat"); library(superheat)}
if(!require(here)){install.packages("here"); library(here)}
if(!require(evolqg)){install.packages("evolqg"); library(evolqg)}
if(!require(doMC)){install.packages("doMC"); library(doMC)}   
if(!require(tictoc)){install.packages("tictoc"); library(tictoc)}   
library(foreach)
registerDoMC(8)

covariates = read_delim(here::here("data/GXEpaper/Covariates_forGEMMA_Jan82021.txt"), 
                        delim = ",")
ID_C = filter(covariates, treatment == 1)$ID
ID_HS = filter(covariates, treatment == 2)$ID

gene_expr = read_delim(here::here("data/GXEpaper/GeneCounts/VOOMCounts_CPM1_head_hsctrl_covfree_4svs_CORRECT_Jan8.21.txt"),
                       delim = "\t", n_max = 100)

gene_expr_list = list(C = select(gene_expr, Gene, all_of(ID_C)), 
                    HS = select(gene_expr, Gene, all_of(ID_HS)))

p = nrow(gene_expr)
x = matrix(NA, p, p)
included = numeric(p)
prop_sig = 0.1
p_cutoff = 1e-4
checkCorrelations = function(i){
    sig = 0
    if(sig <= prop_sig*p){
        for(j in 1:p){
            cr_C = cor.test(as.numeric(gene_expr_list[["C"]][i,-1]), as.numeric(gene_expr_list[["C"]][j,-1]))
            cr_HS = cor.test(as.numeric(gene_expr_list[["HS"]][i,-1]), as.numeric(gene_expr_list[["HS"]][j,-1]))
            p_value = min(cr_C$p.value, cr_HS$p.value)
            if(p_value < p_cutoff) { 
                sig = sig + 1
            }
            if(sig > prop_sig*p) break
        }      
    }
} 
    
tic()
foreach(i = 1:(p-1)) %dopar% {
}
included_genes = included > prop_sig * p
toc()

corr_mat = llply(gene_expr_list, function(x) abs(cor(t(x[,-1]), method = "spearman")))
s_mat = llply(gene_expr_list, function(x) adjacency(t(x[,-1]), power = 6))

#library(patchwork)
#superheat(abs(corr_mat[[1]]), row.dendrogram = T, col.dendrogram = T) 
#superheat(s_mat[[1]], row.dendrogram = T, col.dendrogram = T)

lt = function(x) x[lower.tri(x)]
corr_df = tibble(corr_C = lt(corr_mat[[1]]), 
                 corr_HS = lt(corr_mat[[2]]),
                 s_C = lt(s_mat[[1]]),
                 s_HS = lt(s_mat[[2]]))
library(patchwork)
cc = ggplot(corr_df, aes(corr_C, corr_HS)) + geom_hex(bins = 100) + geom_abline(intercept = 0, slope = 1)
ss = ggplot(corr_df, aes(s_C, s_HS)) + geom_hex(bins = 100) + geom_abline(intercept = 0, slope = 1)
sc = ggplot(corr_df, aes(corr_C, s_C)) + geom_hex(bins = 100)+ geom_abline(intercept = 0, slope = 1)
sc_HS = ggplot(corr_df, aes(corr_HS, s_HS)) + geom_hex(bins = 100)+ geom_abline(intercept = 0, slope = 1)
(cc + ss) / (sc + sc_HS)

cor = corr_mat[[1]]
d = 1 - abs(cor)
a = exp(-d/2)
library(igraph)
g = graph.adjacency(a, weighted = TRUE, mode = 'undirected')
comm = fastgreedy.community(g)
modules = unique(comm$membership)
plot(g)

V(g)$color <- comm$membership+1
g <- set_graph_attr(g, "layout", layout_with_kk(g))
plot(g, vertex.label.dist=1.5)


pickSoftThreshold(t(gene_expr_list[[1]][,-1]) )
pickSoftThreshold(t(gene_expr_list[[2]][,-1]) )

pickHardThreshold(t(gene_expr_list[[2]][,-1]) )


datExprC = t(gene_expr_list[[1]][,-1])
datExprHS = t(gene_expr_list[[2]][,-1])

gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK
sampleTree = hclust(dist(datExpr), method = "average");
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExprHS, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

netC = blockwiseModules(datExprC, power = 2,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       maxBlockSize = 10000,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "C",
                       verbose = 3)
netHS = blockwiseModules(datExprHS, power = 2,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       maxBlockSize = 10000,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "HS",
                       verbose = 3)

sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
