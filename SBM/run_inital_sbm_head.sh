#!/bin/bash

module load conda/4.6.14-1
conda activate gt

python fit_sbm.py --tissue head --layer --wait 10 --graph ../data/output/SBM/graphs/head_fdrLevel-1e-07_genes-3654_density-0.021.xml.gz --output fdr-1e-07 --type nested
python fit_sbm.py --tissue head --layer --wait 10 --graph ../data/output/SBM/graphs/head_fdrLevel-1e-06_genes-4256_density-0.021.xml.gz --output fdr-1e-06 --type nested
python fit_sbm.py --tissue head --layer --wait 10 --graph ../data/output/SBM/graphs/head_fdrLevel-1e-05_genes-4893_density-0.024.xml.gz --output fdr-1e-05 --type nested
python fit_sbm.py --tissue head --layer --wait 10 --graph ../data/output/SBM/graphs/head_fdrLevel-0.0001_genes-5397_density-0.031.xml.gz --output fdr-1e-04 --type nested
python fit_sbm.py --tissue head --layer --wait 10 --graph ../data/output/SBM/graphs/head_fdrLevel-0.001_genes-5567_density-0.049.xml.gz --output fdr-1e-03 --type nested
