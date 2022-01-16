#!/bin/bash

module load conda/4.6.14-1
conda activate gt

python fit_sbm.py --tissue body \
 --graph ../data/output/SBM/graphs/body_fdrLevel-0.001_genes-5124_density-0.059.xml.gz \
  --block ../data/output/SBM/clustering/body_weights-spearman_fdr-1e-03_mcmc_mode_hierarchical-SBM.dill \
 --output fdr-1e-03 \
 --type nested \
 --annealing 100 --wait 100 \
 --mcmc --mcmc-iter 10000
