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

from multipy.fdr import lsu

def filterByEdge(g, corr, cutOff, keepOnlyMain):
    # Filtering edges
    corr = g.edge_properties[corr]
    sign = g.new_ep("bool", True)
    sign.a = np.array(np.abs(corr.a) > cutOff)

    tv = GraphView(g, efilt=sign)

    # Keeping largest component
    if keepOnlyMain:
        comp, hist = label_components(tv)
        main_component = tv.new_vp("bool", (comp.a == np.where(hist == max(hist))[0][0]))
        tv.vertex_properties["main_component"] = main_component
        tv.set_vertex_filter(main_component)
    return tv

def filterByFDR(g, level, keepOnlyMain):
    # Filtering edges
    pvals = np.array(g.edge_properties["pvalue"].a)

    fdr_ep = g.new_ep("bool", True)
    fdr_ep.a = lsu(pvals, q=level)

    tv = GraphView(g, efilt=fdr_ep)

    # Keeping largest component
    if keepOnlyMain:
        comp, hist = label_components(tv)
        main_component = tv.new_vp("bool", (comp.a == np.where(hist == max(hist))[0][0]))
        tv.vertex_properties["main_component"] = main_component
        tv.set_vertex_filter(main_component)
    return tv

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verbose', '-v', action='store_true',
            help='Show more information on the screen.')
    parser.add_argument('--correlation',
            choices=('pearson', 'precision', 'spearman', 'correlation'),
            default='spearman',
            help=("Compute correlation coefficients using either "
                  "'pearson' (standard correlation coefficient), "
                  "'correlation' (shrinkage correlation estimate), "
                  "'precision' (shrinkage inverse correlation estimate), or "
                  "'spearman' (Spearman rank correlation)."))
    parser.add_argument('--fdr', type=float, nargs="*", 
            help=('fdr level used to remove non-significant correlations.'))
    parser.add_argument('--cutoff', type=float,
            help=('The correlation cut-off to use.'))
    parser.add_argument('--tissue', required=True,
            choices=('head', 'body'),
            help='Tissue being analysed.')
    parser.add_argument('--output', '-o', 
            help=('Output file label.'))
    args = parser.parse_args()

    if args.cutoff is None and args.fdr is None:
        raise NameError('Correlation cutoff or fdr must be passed.')
    else:
        print("Reading full graph...")

    if args.tissue == 'body':
        g = load_graph("../data/output/SBM/graphs/VOOMCounts_CPM5_counts4M_covfree_body_ctrl_onlygenesinmainchr_Jul20.21.xml.gz")
    elif args.tissue == 'head':
        g = load_graph("../data/output/SBM/graphs/VOOMCounts_CPM5_counts4M_covfree_head_ctrl_onlygenesinmainchr_Jul20.21.xml.gz")

    # Prune edges by correlation p-value using FDR
    if args.fdr is not None:
        gs = []
        for fdr in args.fdr:
            print("Filtering full graph by fdr level " + str(fdr) + "...")
            tv = filterByFDR(g, fdr, True)
            gs.append(Graph(tv, prune = True))

    # Prune edges by treashold correlation
    if args.cutoff is not None:
        for gi in gs:
            print("Filtering full graph by correlation threshold...")
            tv = filterByEdge(gi, args.correlation, args.cutoff, True)
            gi = Graph(tv, prune = True)

    Ns = []
    densitys = []
    for i, gi in enumerate(gs):
        N = len(gi.get_vertices())
        print(str(N) + " genes")
        Ns.append(N)
        Et = (N * N - N)/2
        E = len(gi.get_edges())
        density = E/Et
        print("Density " + str(density))
        densitys.append(density)

        out_file = "../data/output/SBM/graphs/" + args.tissue
        if args.output is None:
            if args.cutoff is not None:
                out_file = out_file + "_cutoff-" + args.correlation + "_val-" + str(args.cutoff)
            if args.fdr is not None:
                out_file = out_file + "_fdrLevel-" + str(args.fdr[i])
        else:
            out_file = out_file + args.output

        out_file = out_file + '_genes-' + str(N) + '_density-' + str(np.round(density, 3))
        print("Saving filtered graph...")
        print("Output file: " + out_file)
        gi.save(out_file + ".xml.gz")  
        print(" *** ")
    print("Done!")
