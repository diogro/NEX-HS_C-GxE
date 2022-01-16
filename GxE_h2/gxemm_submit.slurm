#!/bin/bash

#sbatch --job-name=gxemm
#sbatch --mem=4G --cpus-per-task=1
#sbatch --time=1:30:00
#sbatch --qos=1wk
#SBATCH --array=1-8786
#SBATCH --output=out/gxemm_%A_%a.out
#SBATCH --error=out/gxemm_%A_%a.err

module load R/4.0.3
cd /Genomics/ayroleslab2/diogro/projects/NEX-HS_C-GxE/R
Rscript GxEMM-HS_C.R --gene $SLURM_ARRAY_TASK_ID --model iid --gene_offset 7999
