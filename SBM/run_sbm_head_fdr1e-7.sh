#!/bin/bash

module load conda/4.6.14-1
conda activate gt

python fit_sbm.py --tissue head --no-layer --wait 10 \
 --graph ../data/output/SBM/graphs/head_fdrLevel-1e-07_genes-1855_density-0.028.xml.gz \
 --output fdr-1e-07 \
 --type nested \
 --block ../data/output/SBM/clustering/head_weights-spearman_fdr-1e-07_mcmc_mode_hierarchical-SBM.dill \
 --mcmc --mcmc-iter 10 
