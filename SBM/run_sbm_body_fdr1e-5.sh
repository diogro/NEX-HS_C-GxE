#!/bin/bash

module load conda/4.6.14-1
conda activate gt

python fit_sbm.py --tissue body \
 --graph ../data/output/SBM/graphs/body_fdrLevel-1e-05_genes-3253_density-0.052.xml.gz \
 --output fdr-1e-05 \
 --type nested \
 --absolute \
 --annealing 100 --wait 100 \
 --mcmc --mcmc-iter 1000 
