#!/bin/bash

module load conda/4.6.14-1
conda activate gt

python fit_sbm.py --tissue body --layer --wait 10 --graph ../data/output/SBM/graphs/body_fdrLevel-1e-07_genes-2212_density-0.046.xml.gz --output fdr-1e-07 --type nested
python fit_sbm.py --tissue body --layer --wait 10 --graph ../data/output/SBM/graphs/body_fdrLevel-1e-06_genes-2717_density-0.044.xml.gz --output fdr-1e-06 --type nested
python fit_sbm.py --tissue body --layer --wait 10 --graph ../data/output/SBM/graphs/body_fdrLevel-1e-05_genes-3407_density-0.042.xml.gz --output fdr-1e-05 --type nested
python fit_sbm.py --tissue body --layer --wait 10 --graph ../data/output/SBM/graphs/body_fdrLevel-0.0001_genes-4432_density-0.039.xml.gz --output fdr-1e-04 --type nested
python fit_sbm.py --tissue body --layer --wait 10 --graph ../data/output/SBM/graphs/body_fdrLevel-0.001_genes-5443_density-0.043.xml.gz --output fdr-1e-03 --type nested
