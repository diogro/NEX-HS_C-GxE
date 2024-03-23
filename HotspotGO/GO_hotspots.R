source(here::here("R/go_functions.R"))
library(furrr)

plots_path = "B:/Dropbox/labbio/articles/NEX_BodyHead_Control-SBM/figures/"

# Read the genes.csv file and extract the first column as the background for the GO analysis
genes_background <-
  list(
    body = read_csv(here::here("HotspotGO/bodygenebackground.txt"), 
                     col_names = "Genes") |> 
            pull(Genes) |> 
            unique(),
    head = read_csv(here::here("HotspotGO/headgenebackground.txt"),
                     col_names = "Genes") |>
            pull(Genes) |> 
            unique()
    )

data_head = read_delim(here::here("HotspotGO/HEADhotspots_targetgenes_forGO_6Dec23.txt"))
data_body = read_delim(here::here("HotspotGO/BODYhotspots_targetgenes_forGO_6Dec23.txt"))

# genes = data_head[["X3L_20369028"]]
tissue = "head"
perform_enrichment(data_head[["X3L_20369028"]], "head")
perform_enrichment <- function(genes, tissue) {
  genes = na.omit(genes)
  go_enrichment_result <- enrichGO(
    gene = genes,
    universe = genes_background[[tissue]],
    OrgDb = org.Dm.eg.db,
    keyType = "FLYBASE",
    ont = "BP",  # Biological Process ontology
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
  )
  
 return(go_enrichment_result)
}

# Use purrr's map function to perform enrichment for each columns
plan(multisession, workers = 16)
enResults_body <- names(data_body) |> 
  set_names() |>
  future_map(\(x) perform_enrichment(pull(data_body, x), 
                              "body"))

enResults_head <- names(data_head) |> 
  set_names() |>
  map(\(x) perform_enrichment(pull(data_head, x), 
                              "head"))

enTables = list(body = ldply(enResults_body, 
                      function(x) x@result) %>%
  filter(p.adjust < 0.05, Count >= 4),
head = ldply(enResults_head, 
                      function(x) x@result) %>%
  filter(p.adjust < 0.05, Count >= 4)) |> 
    list_rbind(names_to = "Tissue") |>
    rename(Position = .id)


# Write the results to a file
write_csv(enTables, pointblank::affix_date(here::here("HotspotGO/HotSpotGoEnrichmentTable.csv")))

enriched_pos = enTables |>
  group_by(Tissue, Position) |>
  count(name = "terms")

all_pos = c(names(data_body), names(data_head)) 
length(which(all_pos %in% enriched_pos$Position)) / length(all_pos)
