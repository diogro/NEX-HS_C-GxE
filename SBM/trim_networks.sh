#!/bin/bash

module load conda/4.6.14-1
conda activate gt

python ./trim_networks.py --tissue head --fdr 0.05 0.01 0.001 0.0001 0.00001 0.000001 0.0000001 
python ./trim_networks.py --tissue body --fdr 0.05 0.01 0.001 0.0001 0.00001 0.000001 0.0000001 