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

from graph_tool.all import *
import pandas as pd
import numpy as np
import scipy as sp
from sklearn.covariance import LedoitWolf, OAS
import statsmodels.api as sm

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verbose', '-v', action='store_true',
            help='Show more information on the screen.')
    parser.add_argument('--data', required=True,
            help=('Path to the covfree expression data.'))
    parser.add_argument('--out', required=True,
            help=('Path to the output graph file.'))
    
    args = parser.parse_args()
    
    if args.verbose:
        print('sys.argv:')
        print(sys.argv)
        print()
        print('numpy version:', np.__version__)
        print('pandas version:', pd.__version__)
        print('scipy version:', sp.__version__)
        print()
    
    gene_expr_raw = pd.read_table(args.data)
    gene_expr = gene_expr_raw.T

    X_centered = (gene_expr - gene_expr.mean()) / np.sqrt(gene_expr.var())
    oa = OAS(store_precision=True, assume_centered=True)
    gene_expr_OAS_corr = oa.fit(X_centered)

    n_genes = gene_expr_OAS_corr.covariance_.shape[1]
    g = Graph(directed=False)
    g.add_vertex(n = n_genes)

    spearman = g.new_ep("double", 0)
    pval = g.new_ep("double", 0)
    genes = g.new_vertex_property("string", np.array(np.array(gene_expr.columns, dtype = "str")))
    g.vertex_properties["genes"] = genes

    for i in range(n_genes):
        for j in range(i):
            spearman_r = sp.stats.spearmanr(X_centered.iloc[:,i], X_centered.iloc[:,j])
            g.add_edge(i, j)
            e = g.edge(i, j)
            pval[e] = spearman_r[1]
            spearman[e] = spearman_r[0]

    g.edge_properties["pvalue"] = pval
    g.edge_properties["spearman"] = spearman

    g.save(args.out)