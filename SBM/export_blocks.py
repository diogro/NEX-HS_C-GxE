import os
os.environ["OMP_NUM_THREADS"] = "32"

from contextlib import contextmanager
import argparse
import os.path
import csv
import time
import sys
from functools import partial
import shutil as sh

import dill

from graph_tool.all import *
import pandas as pd
import numpy as np
import scipy as sp
from sklearn.covariance import LedoitWolf, OAS
import statsmodels.api as sm

from multipy.fdr import lsu

def ensure_dir(file_path):
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)
    if not os.path.exists(file_path):
        os.makedirs(file_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verbose', '-v', action='store_true',
            help='Show more information on the screen.')
    parser.add_argument('--tissue', required=True,
            choices=('head', 'body'),
            help='Tissue being analysed.')
    parser.add_argument('--label', required=True, 
            help=('Nested block partition csv and dill output file.'))
    parser.add_argument('--path', dest='input_path', default = '../data/output/SBM/clustering/',
            help=('Nested block partition csv and dill output file.'))
    parser.add_argument('--out', dest='out_path', default = None,
            help=('Outputh folder name.'))
    args = parser.parse_args()

    block_df = pd.read_csv(args.input_path + args.label + ".csv")

    with open (args.input_path + args.label + ".dill", "rb") as fh:
        emp = dill.load(fh)
    g = emp[0]
    bs = emp[1]
    print("Read dilled blockmodel...")
    state = minimize_nested_blockmodel_dl(g, init_bs=bs, 
                                          state_args=dict(recs=[g.ep.z_s],
                                                          rec_types=["real-normal"]))
    print("Recreated blockmodel...")

    if args.out_path is None:
        out_folder = args.input_path + args.label + '_gene-blocks'
    else:
        out_folder = args.out_path
    ensure_dir(out_folder)

    print("Clearing output folder...")
    for filename in os.listdir(out_folder):
        file_path = os.path.join(out_folder, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                sh.rmtree(file_path)
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))

    print("Calculating block sizes...")
    blocks = [list(set(block_df[b])) for b in block_df.filter(like='B', axis=1)]
    block_sizes = [len(b) for b in blocks]
    block_sizes = -np.sort(-np.array(list(set(block_sizes))))
    block_sizes = [x for x in block_sizes if x >= 2]

    block_df["Gene"].to_csv(out_folder + "/background.csv", header=False, index=False )

    print("Creating gene lists...")
    n_levels = len(block_sizes)
    output_df = pd.DataFrame(columns=('Nested_Level', 'Block', 'File', 'N_genes', 'Internal_degree', 'Assortatitvity'))
    l = 0
    for i in range(n_levels):
        print("At level: " + str(i+1))
        bl = blocks[i]
        for b in bl:

            line = [i+1]
            line.append(b)

            df = block_df[block_df['B' + str(i+1)]==b]
            genes = df["Gene"]
            file_name = "/" + '-'.join([str(num) for num in list(df.filter(like='B', axis=1).iloc[1,range(i, n_levels)])]) + ".csv"

            line.append(file_name)
            line.append(genes.shape[0])

            ensure_dir(out_folder + "/Level_" + str(i+1))
            genes.to_csv(out_folder + "/Level_" + str(i+1) + file_name, header=False, index=False )

            # Weighted
            ers = adjacency(state.levels[i].bg, weight=state.levels[i].mrs)
            B = state.levels[i].bg.num_vertices()
            M = ers.sum()
            q_r = B * (ers[b,b] - ers[b,:].sum() ** 2/M)/M 

            line.append(ers[b,b])
            line.append(q_r)

            output_df.loc[l] = line
            l = l + 1
    output_df.to_csv(out_folder + "/block_summary.csv", index=False )


    