library(tidyverse)

label = "head_weights-spearman_fdr-1e-02"

addRawDegree <- function(label){
    gene_level = read_csv(paste0("/Genomics/argo/users/damelo/projects/NEX-HS_C-GxE/data/output/SBM/clustering/", label, "_mcmc_mode_hierarchical-SBM.csv"))

    block_summary_file = paste0("/Genomics/argo/users/damelo/projects/NEX-HS_C-GxE/data/output/SBM/clustering/", label, "_mcmc_mode_hierarchical-SBM_gene-blocks/block_summary.csv")
    block_level = read_csv(block_summary_file)

    average_degrees = list()
    for(level in 1:5){
        B = paste0("B", level)
        average_degrees[[level]] = gene_level %>%
            group_by_({{B}}) %>%
            summarise(Average_Spearman_degree = mean(E_corr)) %>%
            rename(Block = {{B}}) %>%
            mutate(Nested_Level = level)
    }
    inner_join(block_level, do.call(bind_rows, average_degrees), by = c("Block", "Nested_Level"))
}
addRawDegree("body_weights-spearman_fdr-1e-03")
block_level = addRawDegree("head_weights-spearman_fdr-1e-02")
block_level |> as.data.frame() |> head()
as.data.frame(block_level)[1,] 
block_level[,4:10] |> cor()
