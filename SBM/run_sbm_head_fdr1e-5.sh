#!/bin/bash

module load conda/4.6.14-1
conda activate gt

python fit_sbm.py --tissue head --no-layer --wait 2 \
 --graph ../data/output/SBM/graphs/head_fdrLevel-1e-05_genes-2777_density-0.023.xml.gz \
 --output fdr-1e-05 \
 --type nested \
 --block ../data/output/SBM/clustering/head_weights-spearman_fdr-1e-05_mcmc_mode_hierarchical-SBM.dill \
 --mcmc --mcmc-iter 10000
