#!/bin/bash

module load conda/4.6.14-1
conda activate gt

python fit_sbm.py --tissue head \
 --graph ../data/output/SBM/graphs/head_fdrLevel-0.0001_genes-3589_density-0.021.xml.gz \
 --output fdr-1e-04 \
 --type nested \
 --absolute \
 --annealing 100 --wait 100 \
 --mcmc --mcmc-iter 1000
