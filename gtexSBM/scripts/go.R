source(here::here("gtexSBM/scripts/go_functions_gtex.R"))
library(doMC)
registerDoMC(32)
go_set = makeEnrichment(block_path = here::here("data/output/SBM/gtex/blockSummary/fdr-1e-3/COLON"))

table_en = table(go_set$summary$n_enrich[go_set$summary$Nested_Level==4]!=0)
table_en
table_en/sum(table_en)

x = go_set[[1]]
names(x)
