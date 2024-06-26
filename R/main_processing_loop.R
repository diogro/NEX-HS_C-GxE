main_count_processing <- function(dset_name,
                                  exp_data,
                                  experimental_metadata,
                                  feature_vec,
                                  assay_name,
                                  min_expr = 1) {
  tic()              
  print(dset_name)
  dset <- exp_data[[dset_name]]

  # Catch empty columns_to_ignore
  columns_to_ignore = tryCatch(
    expr = {
      unlist(strsplit(experimental_metadata[experimental_metadata$id == dset_name, ]$columns_to_ignore, split = ";"))
    },
    error = function(e) {
      columns_to_ignore <- c("")
      columns_to_ignore
    }
  )
  columns_to_ignore = c(columns_to_ignore, crap_cols)

  metadata <- colData(dset)
  metadata <- data.frame(metadata)

  counts <- assays(dset)[[assay_name]]

  if(dset_name == "SKIN"){
      remove_cell_culture = metadata[["gtex.smtsd"]] != "Cells - Cultured fibroblasts"
      metadata = metadata[ remove_cell_culture,]
      counts   =   counts[, remove_cell_culture]
  }
  if(dset_name == "ESOPHAGUS"){
      remove_cell_culture = metadata[["gtex.smtsd"]] == "Esophagus - Mucosa"
      metadata = metadata[ remove_cell_culture,]
      counts   =   counts[, remove_cell_culture]
  }

  print(paste0("Unfiltered count dimensions: ", dim(counts)[1], " x ", dim(counts)[2]))
  print(paste0("Unfiltered metadata dimensions: ", dim(metadata)[1], " x ", dim(metadata)[2]))

  filtered_data <- make_filtered_data(counts, metadata, feature_vec)
  metadata <- filtered_data$metadata
  counts <- filtered_data$counts
  n_samples <- dim(metadata)[1]
  metadata$sample_id <- rownames(metadata)

  print("Normalizing and estimating mean-variance weights...")
  countdata.list <- DGEList(counts = counts, samples = metadata, genes = rownames(dset))
  rep_names = c("technical_replicate_group", "wells.replicate", "individual")
  if (any(names(metadata) %in% rep_names)) {
    print("Summing technical replicates...")
    rep_col = names(metadata)[names(metadata) %in% rep_names]
    countdata.list <- sumTechReps(countdata.list, metadata[[rep_col]])
  }

  print("Calculating normalization factors...")
  countdata.norm <- calcNormFactors(countdata.list)

  print("Trimming...")
  cutoff <- inv_log2_plus05(min_expr)
  drop <- which(apply(cpm(countdata.norm), 1, max) < cutoff)
  countdata.norm <- countdata.norm[-drop, ]

  drop <- which(apply(cpm(countdata.norm), 1, min) < 1)
  countdata.norm <- countdata.norm[-drop, ]

  drop <- which(apply(cpm(countdata.norm), 1, mean) < 5)
  countdata.norm <- countdata.norm[-drop, ]

  # recalc with updated metadata
  values_count <- sapply(lapply(countdata.norm$samples, unique), length)
  countdata.norm$samples <- countdata.norm$samples[, names(countdata.norm$samples[, values_count > 1])]

  print("Making design matrix...")
  countdata.norm$samples  <- remove_redundant_features(countdata.norm$samples )
  # Removing columns with a crazy number of levels that mess everything up.
  # (this is why we have random effects by the way)
  countdata.norm$samples  <- remove_large_factors(
    countdata.norm$samples,
    columns_to_ignore
  )
  countdata.norm$samples <- select_meta(countdata.norm$samples)
  design <- make_design_matrix(countdata.norm$samples, columns_to_ignore)
  print(paste("Design matrix size:", paste(dim(design), collapse = " x ")))

  pca_on_raw <- pca_plot(countdata.norm$counts, color = rep("1", ncol(countdata.norm$counts)))
  #pca_on_raw <- pca_plot(countdata.norm$counts, color = metadata[["gtex.smtsd"]])
  #save_plot("pca_on_raw.png", pca_on_raw, base_height = 6)
  #screen_on_raw <- scree_plot(countdata.norm$counts)

  if(dset_name == "BLOOD"){
      bigones = sort(apply(countdata.norm$counts, 1, max), decreasing = T)
      remove_genes = which(rownames(countdata.norm) %in% names(bigones)[1:3])
      countdata.norm  =   countdata.norm[-remove_genes,]
  }
  if(dset_name %in% c("COLON", "STOMACH")){
      bigones = sort(apply(countdata.norm$counts, 1, max), decreasing = T)
      remove_genes = which(rownames(countdata.norm) %in% names(bigones)[1])
      countdata.norm  =   countdata.norm[-remove_genes,]
  }
  print(paste0("Filtered count dimensions: ", dim(countdata.norm$counts)[1], " x ", dim(countdata.norm$counts)[2]))
  countdata_resids <- cpm_lm(countdata.norm, design = design)

  countdata.norm_noOut <- countdata.norm
  
  rpca_resid <- PcaGrid(t(countdata_resids), 20, crit.pca.distances = 0.99)
  countdata.norm_noOut$counts <- countdata.norm_noOut$counts[, rpca_resid@flag]
  countdata.norm_noOut$samples <- countdata.norm_noOut$samples[rpca_resid@flag, ,drop = FALSE]

  # PCA plot with Batch effects (this plot happens here to make use of the outlier tags from the robust PCA)
  pca_on_resids <- pca_plot(countdata_resids, color = !rpca_resid@flag)
  #scree_on_resids <- scree_plot(countdata_resids)

  print(paste0(
    "Filtered metadata dimensions: ",
    paste(dim(countdata.norm_noOut$samples), collapse = " x ")
  ))

  design_noOut <- make_design_matrix(countdata.norm_noOut$samples, columns_to_ignore)
  countdata_resids_noOut <- cpm_lm(countdata.norm_noOut, design = design_noOut)
  pca_on_resids_noOut <- pca_plot(countdata_resids_noOut, color = rep("1", ncol(countdata_resids_noOut)))

  # Batch effects With PC1

  PCs <- pca(countdata_resids_noOut)
  countdata.norm_noOut$samples$PC1 <- PCs$pc[1, ]
  design_with_pc1 <- make_design_matrix(countdata.norm_noOut$samples, columns_to_ignore)

  countdata_resids_with_pc1 <- cpm_lm(countdata.norm_noOut, design = design_with_pc1)

  # PCA plot with Batch effects and PC1
  pca_on_resids_with_pc1 <- pca_plot(countdata_resids_with_pc1, color = rep("1", ncol(countdata_resids_with_pc1)))

  # print("Writing figures")

  plt <- plot_grid(
    nrow = 2, scale = 0.9,
    pca_on_raw + ggtitle("Uncorrected"),
    pca_on_resids + ggtitle("Known effects"),
    pca_on_resids_noOut + ggtitle("Known effects, no outliers"),
    pca_on_resids_with_pc1 + ggtitle("Known effects and PC1")
  ) 
  time = toc()
  # print("Appending results and metadata to lists")
  list(
    name = dset_name,
    n_samples = dim(countdata.norm_noOut$samples)[1],
    normcounts_raw = countdata.norm,
    normcounts_noOut = countdata.norm_noOut,
    plotPanel = plt,
    plot_list = list(uncorrected = pca_on_raw, batch = pca_on_resids, clean = pca_on_resids_noOut, pc1 = pca_on_resids_with_pc1),
    residuals_noOut = countdata_resids_noOut,
    residuals_raw = countdata_resids,
    residuals_pc1 = countdata_resids_with_pc1,
    metadata = countdata.norm_noOut$samples,
    time = time
  )
}

main_loop <- function (...) {
  return(tryCatch(main_count_processing(...), error=function(e) NA))
}
