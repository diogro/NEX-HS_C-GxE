#!/bin/bash

module load conda/4.6.14-1
conda activate gt

python fit_sbm.py --tissue head \
 --graph ../data/output/SBM/graphs/head_fdrLevel-0.01_genes-5261_density-0.033.xml.gz \
 --block ../data/output/SBM/clustering/head_weights-spearman_fdr-1e-02_mcmc_mode_hierarchical-SBM.dill \
 --output fdr-1e-02 \
 --type nested \
 --annealing 100 --wait 100 \
 --mcmc --mcmc-iter 10000
