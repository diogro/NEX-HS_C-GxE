#!/bin/bash

module load conda/4.6.14-1
conda activate gt

python fit_sbm.py --tissue body --no-layer --wait 10 --graph ../data/output/SBM/graphs/body_fdrLevel-1e-07_genes-2314_density-0.045.xml.gz --output fdr-1e-07 --type nested
python fit_sbm.py --tissue body --no-layer --wait 10 --graph ../data/output/SBM/graphs/body_fdrLevel-1e-06_genes-2750_density-0.047.xml.gz --output fdr-1e-06 --type nested
python fit_sbm.py --tissue body --no-layer --wait 10 --graph ../data/output/SBM/graphs/body_fdrLevel-1e-05_genes-3253_density-0.052.xml.gz --output fdr-1e-05 --type nested
python fit_sbm.py --tissue body --no-layer --wait 10 --graph ../data/output/SBM/graphs/body_fdrLevel-0.0001_genes-4033_density-0.055.xml.gz --output fdr-1e-04 --type nested
python fit_sbm.py --tissue body --no-layer --wait 10 --graph ../data/output/SBM/graphs/body_fdrLevel-0.001_genes-5124_density-0.059.xml.gz --output fdr-1e-03 --type nested
