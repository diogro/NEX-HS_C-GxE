#!/bin/bash

module load conda/4.6.14-1
conda activate gt

python fit_sbm.py --tissue body --no-layer --wait 10 \
 --graph ../data/output/SBM/graphs/body_fdrLevel-1e-06_genes-2750_density-0.047.xml.gz \
 --output fdr-1e-06 \
 --type nested \
 --block ../data/output/SBM/clustering/body_weights-spearman_fdr-1e-06_hierarchical-SBM.dill \
 --mcmc --mcmc-iter 1000 
