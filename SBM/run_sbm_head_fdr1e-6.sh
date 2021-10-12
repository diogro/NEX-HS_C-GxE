#!/bin/bash

module load conda/4.6.14-1
conda activate gt

python fit_sbm.py --tissue head --no-layer --wait 10 \
 --graph ../data/output/SBM/graphs/head_fdrLevel-1e-06_genes-2234_density-0.026.xml.gz \
 --output fdr-1e-06 \
 --type nested \
 --block ../data/output/SBM/clustering/head_weights-spearman_fdr-1e-06_hierarchical-SBM.dill \
 --mcmc --mcmc-iter 1000 
