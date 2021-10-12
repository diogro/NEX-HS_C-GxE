#!/bin/bash

module load conda/4.6.14-1
conda activate gt

python fit_sbm.py --tissue head --no-layer --wait 100 --graph ../data/output/SBM/graphs/head_fdrLevel-1e-07_genes-1855_density-0.028.xml.gz --output fdr-1e-07 --type nested
python fit_sbm.py --tissue head --no-layer --wait 100 --graph ../data/output/SBM/graphs/head_fdrLevel-1e-06_genes-2234_density-0.026.xml.gz --output fdr-1e-06 --type nested
python fit_sbm.py --tissue head --no-layer --wait 100 --graph ../data/output/SBM/graphs/head_fdrLevel-1e-05_genes-2777_density-0.023.xml.gz --output fdr-1e-05 --type nested
python fit_sbm.py --tissue head --no-layer --wait 100 --graph ../data/output/SBM/graphs/head_fdrLevel-0.0001_genes-3589_density-0.021.xml.gz --output fdr-1e-04 --type nested
python fit_sbm.py --tissue head --no-layer --wait 100 --graph ../data/output/SBM/graphs/head_fdrLevel-0.001_genes-4685_density-0.02.xml.gz --output fdr-1e-03 --type nested
