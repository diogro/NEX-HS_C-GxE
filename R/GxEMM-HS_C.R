library(optparse)
library(bracer)

if (sys.nframe() == 0L) {
  option_list <- list(
    make_option("--gene",
                help = ("Gene number to run"),
                metavar = "i_gene"),
    make_option("--model",
                help = ("Model to run in GxEMM"),
                default = "full",
                metavar = "model"),
    make_option("--input",
                help = ("Path to gene expression file"),
                default = "/Genomics/ayroleslab2/lamaya/bigProject/Datasets/GXEpaper/GeneCounts/VOOMCounts_CPM1_head_hsctrl_covfree_4svs_CORRECT_Jan8.21.txt",
                metavar = "input_path"),
    make_option("--grm", default = "/Genomics/ayroleslab2/lamaya/bigProject/Datasets/GXEpaper/GRM/head.hsc.final.norelatedness1.miss50.sXX.txt",
                help = ("Path to GRM file"),
                metavar = "GRM_path"),
    make_option("--covariates", default = "/Genomics/ayroleslab2/lamaya//bigProject/Datasets/GXEpaper/Covariates_forGEMMA_Jan82021.txt",
                help = ("Path to covariates file"),
                metavar = "covariates_path"),
    make_option("--outputDir",
                default = "/Genomics/ayroleslab2/diogro/projects/NEX-HS_C-GxE/data/output/GxEMM/",
                help = ("OutPut directory"),
                metavar = "outputDir")
  )
  parser_object <- OptionParser(usage = "Rscript %prog [Options] --gene [gene]\n",
                                option_list = option_list,
                                description = "Run GxEMM")
  
  ## TO TEST INTERACTIVELY the command-line arguments
  #input <- "--dir ../dados/estado_SP/SRAG_hospitalizados/dados/ --escala municipio --geocode 350750 --dataBase 2020_05_20"
  #command.args <- strsplit(input, " ")[[1]]
  #opt <- parse_args(parser_object, args = command.args, positional_arguments = TRUE)
  ## SKIP opt line below
  ## aliases
  opt <- parse_args(parser_object, args = commandArgs(trailingOnly = TRUE),
                    positional_arguments = TRUE)
  i_gene <- opt$options$i_gene
  model <- opt$options$model
  input_path <- opt$options$input_path
  GRM_path <- opt$options$GRM_path
  covariates_path <- opt$options$covariates_path
  out.dir <- opt$options$outputDir
  
  # quit on error when run non-interactively #small change because this is killing my local sessions T_T
  if (!interactive()) options(error = function() quit(save = "no", status = 1))
}

library(GxEMM)
library(stringr)
ldak_loc  <- "/Genomics/grid/users/damelo/.local/bin/ldak5.linux "

expr_df = t(read.table(input_path, row.names = 1))
rownames(expr_df) = gsub("^X", "", rownames(expr_df))

covariates = read.table(covariates_path)
ID_C = row.names(covariates[covariates$treatment == 1,])
ID_HS = row.names(covariates[covariates$treatment == 2,])

Z_HS = as.numeric(rownames(expr_df) %in% ID_HS)
Z_C = as.numeric(rownames(expr_df) %in% ID_C)

Z_all = cbind(Z_HS, Z_C)
X     <- Z_all[,-1] 

GRM = as.matrix(read.table(GRM_path))

out_free_HSC = GxEMM(expr_df[,i_gene], X = X, K = GRM, Z = Z_all, 
                     gtype='free', etype='free', ldak_loc=ldak_loc)
out_file = file.path(out.dir, paste0(str_pad(i_gene, 4, pad = "0"), "_", colnames(expr_df)[i_gene], "_", model, ".rds"))
saveRDS(out_free_HSC, file = out_file)
