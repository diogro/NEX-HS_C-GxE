#!/bin/bash

module load conda/4.6.14-1
conda activate gt

python fit_sbm.py --tissue head --no-layer --wait 0 \
 --graph ../data/output/SBM/graphs/head_fdrLevel-0.0001_genes-3589_density-0.021.xml.gz \
 --output fdr-1e-04 \
 --type nested \
 --block ../data/output/SBM/clustering/head_weights-spearman_fdr-1e-04_mcmc_mode_hierarchical-SBM.dill \
 --mcmc --mcmc-iter 10000
