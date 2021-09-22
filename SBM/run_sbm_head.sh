#!/bin/bash

module load conda/4.6.14-1
conda activate gt

values=(0.55 0.5 0.45 0.40 0.35 0.3)
for i in "${values[@]}"
    do
        python fit_sbm.py \
            --data VOOMCounts_CPM5_counts4M_covfree_head_ctrl_onlygenesinmainchr_Jul20.21_regularized_correlations_precisions_spearman_correlation_cutoff_0.1.xml.gz \
            --correlation spearman \
            --tissue head \
            --sigma $i \
            --type all \
            --wait 1000 
    done