
source("~/projects/exp_var/functions.R")

library(tictoc)

# Set timeout to avoid failure when trying to download
# GTEx or other large datasets
options(timeout = 1800)

# Load metadata from recount3
experimental_metadata_rc3 <- read_csv("~/projects/exp_var/recount3_metadata.csv")[1:17,]

plots_dir <- here::here("R/plots/")

# Setup filters for removing samples accordingly to the metadata
feature_vec <- list()
feature_vec[["disease"]] <- c("normal", "control", "", NA,
                              "non inflammatory bowel disease control")
feature_vec[["treatment"]] <- c("normal", "control", "", NA)

# This is commented out as we're sticking with Expression Atlas for now.
library(recount3)

# Move cache to deal with quota issues
cache <- recount3_cache(cache_dir = here::here("R/cache"))
human_projects <- available_projects(bfc = cache)

source(here::here("./R/main_processing_loop.R"))
parallel <- FALSE
if (.Platform$OS.type == "unix") {
  parallel <- TRUE
  library(doMC)
  registerDoMC(64)
}

pull_data <- TRUE
if (pull_data) {
tic()
  gtex_list = llply(experimental_metadata_rc3$id, downloadRecount3, .parallel = parallel)
  names(gtex_list) = experimental_metadata_rc3$id
toc()
}

results_list_gtex <- llply(names(gtex_list),
                      main_count_processing,
                      exp_data = gtex_list,
                      experimental_metadata = experimental_metadata_rc3,
                      feature_vec = feature_vec,
                      assay_name = "raw_counts",
                      min_expr = 1,
                      .parallel = parallel)
names(results_list_gtex) <- names(gtex_list)
save(results_list_gtex, file = here::here("R/cache/results_list_gtex.RData"))

for (dset_name in names(results_list_gtex)){
  save_plot(filename = paste0(plots_dir, dset_name, "_pca.png"),
            results_list_gtex[[dset_name]]$plotPanel,
            base_height = 6, base_asp = 1.2, ncol = 2, nrow = 2)
}

csv_path = here::here("data/output/SBM/gtex/raw")
sbm_output = vector("list", length(results_list_gtex))
names(sbm_output) = names(results_list_gtex)
for(i in 1:length(sbm_output)){
  x = results_list_gtex[[i]]
  sbm_output[[i]] = as_tibble(x$residuals_pc1) %>% 
    dplyr::mutate(Gene = rownames(x$residuals_pc1)) %>% 
    dplyr::select(Gene, everything())
}
laply(sbm_output, dim)
gene_df = Reduce(function(x, y) inner_join(dplyr::select(x, Gene), dplyr::select(y, Gene)),
  sbm_output)
sbm_output = llply(sbm_output, inner_join, gene_df)
laply(sbm_output, dim)
for(i in 1:length(sbm_output)){
  name = names(sbm_output)[i]
  x = as.data.frame(sbm_output[[i]])
  rownames(x) = x$Gene
  x$Gene = NULL
  write.table(x, file = file.path(csv_path, paste0(i, "_", name, ".csv")), sep = "\t")
  print(name)
}
