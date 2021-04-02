if(!require(WGCNA)){BiocManager::install("WGCNA"); library(WGCNA)}
if(!require(plyr)){install.packages("plyr"); library(plyr)}
if(!require(tidyverse)){install.packages("tidyverse"); library(tidyverse)}
if(!require(superheat)){install.packages("superheat"); library(superheat)}

covariates = read_delim("./data/GXEpaper/Covariates_forGEMMA_Jan82021.txt", 
                        delim = ",")
ID_C = filter(covariates, treatment == 1)$ID
ID_HS = filter(covariates, treatment == 2)$ID

gene_expr = read_delim("./data/GXEpaper/GeneCounts/VOOMCounts_CPM1_head_hsctrl_covfree_4svs_CORRECT_Jan8.21.txt", delim = "\t", 
                       n_max = 500)

gen_expr_list = list(C = select(gene_expr, Gene, all_of(ID_C)), 
                    HS = select(gene_expr, Gene, all_of(ID_HS)))

corr_mat = llply(gen_expr_list, function(x) cor(t(x[,-1])))
superheat(corr_mat[[1]], row.dendrogram = T, col.dendrogram = T)
superheat(corr_mat[[2]], row.dendrogram = T, col.dendrogram = T)


