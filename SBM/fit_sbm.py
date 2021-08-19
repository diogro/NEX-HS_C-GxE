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

def filterByEdge(g, corr, cutOff, keepOnlyMain):
    # Filtering edges
    corr = g.edge_properties[corr]
    sign = g.new_ep("bool", True)
    sign.a = np.array(corr.a > cutOff)

    tv = GraphView(g, efilt=sign)

    # Keeping largest component
    if keepOnlyMain:
        comp, hist = label_components(tv)
        main_component = tv.new_vp("bool", (comp.a == np.where(hist == max(hist))[0][0]))
        tv.vertex_properties["main_component"] = main_component
        tv.set_vertex_filter(main_component)
    return tv

def create_Block_df(g, corr, state):
    genes = g.vertex_properties["genes"]
    block_df = pd.DataFrame(columns=('Gene', "Degree", "E_corr", 'Block'))
    for v in g.vertex_index:
        line = [genes[v]]
        line.append(g.get_total_degrees([v])[0])
        line.append(np.mean(g.get_all_edges(v, [corr] )[:,2]))
        line.append(state.get_blocks()[v])
        block_df.loc[v] = line
    return block_df

def get_group(x, state):
    levels = state.get_levels()
    n_levels = 5#len(levels)
    r = np.zeros(n_levels)
    r[0] = levels[0].get_blocks()[x]
    for i in range(1, n_levels):
        r[i] = levels[i].get_blocks()[r[i-1]]
    r = r.astype(int)
    return r
def create_nestedBlock_df(g, corr, state):
    genes = g.vertex_properties["genes"]
    nested_block_df = pd.DataFrame(columns=('Gene', "Degree", "E_corr", 'B1', "B2", "B3", "B4", "B5"))
    for v in g.vertex_index:
        line = [genes[v]]
        line.append(g.get_total_degrees([v])[0])
        line.append(np.mean(g.get_all_edges(v, [corr] )[:,2]))
        [line.append(i) for i in get_group(v, state)]
        nested_block_df.loc[v] = line
    return nested_block_df

def run_non_nested_SBM(g, corr, args):
    n_trials = args.trials
    state_nn_min_list = []
    print("Running SBM " + str(n_trials) + " times")
    for i in range(n_trials):
        print(i)
        #state_min = minimize_blockmodel_dl(gf, state_args=dict(recs=[corr],
        #                                                       rec_types=["real-normal"]))
        state_min = minimize_blockmodel_dl(g)
        state_nn_min_list.append(state_min)

    print("Improve solution with merge-split")
    ret_nn = []
    state_nn_mcmc_list = []
    for j in range(n_trials):
        print(j)
        for i in range(200):
            state = state_nn_min_list[j].copy(sampling = True)
            x = state.multiflip_mcmc_sweep(niter=20, beta=np.inf)
            state_nn_mcmc_list.append(state)
            ret_nn.append(x)

    description_lenghts = np.zeros(n_trials)
    for j in range(n_trials):
        description_lenghts[j] = state_nn_mcmc_list[j].entropy()
    min_index = np.argmin(description_lenghts)
    state_min_nn = state_nn_mcmc_list[min_index]

    plot_file = "../data/output/SBM/plots/graph_plot_" + args.tissue + "_cutoff-" + args.correlation + "_val-" + str(args.sigma) + "_clustered-non-hierarchical-SBM.png"
    state_min_nn.draw(output = plot_file)

    block_df = create_Block_df(g, corr, state_min_nn)
    output_file = "../data/output/SBM/clustering/" + args.tissue + "_cutoff-" + args.correlation + "_val-" + str(args.sigma) + "_non-hierarchical-SBM.csv"
    block_df.to_csv(output_file)

def run_nested_SBM(g, corr, args):
    n_trials = args.trials
    print("Running Hierarchical SBM " + str(n_trials) + " times")
    state_min_list = []
    for i in range(n_trials):
        print(i)
        #state_min = minimize_nested_blockmodel_dl(gf, state_args=dict(recs=[corr],
        #                                                              rec_types=["real-normal"]))
        state_min = minimize_nested_blockmodel_dl(g)
        state_min_list.append(state_min)

    # improve solution with merge-split
    # state_min = state_min.copy(bs=state_min.get_bs() + [np.zeros(1)] * 4, sampling=True)

    ret = []
    state_mcmc_list = []
    print("Improve solution with merge-split")
    for j in range(n_trials):
        print(j)
        for i in range(200):
            state = state_min_list[j].copy(sampling = True)
            x = state.multiflip_mcmc_sweep(niter=20, beta=np.inf)
            state_mcmc_list.append(state)
            ret.append(x)

    description_lenghts = np.zeros(n_trials)
    for j in range(n_trials):
        description_lenghts[j] = state_mcmc_list[j].entropy()
    min_index = np.argmin(description_lenghts)
    state_min = state_mcmc_list[min_index]

    levels = state_min.get_levels()

    plot_file = "../data/output/SBM/plots/graph_plot_" + args.tissue + "_cutoff-" + args.correlation + "_val-" + str(args.sigma) + "_clustered-hierarchical-SBM.png"
    state_min.draw(output = plot_file)

    nested_block_df = create_nestedBlock_df(g, corr, state_min)
    output_file = "../data/output/SBM/clustering/" + args.tissue + "_cutoff-" + args.correlation + "_val-" + str(args.sigma) + "_hierarchical-SBM.csv"
    nested_block_df.to_csv(output_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verbose', '-v', action='store_true',
            help='Show more information on the screen.')
    parser.add_argument('--correlation',
            choices=('pearson', 'precision', 'spearman', 'correlation'),
            default='pearson',
            help=("Compute correlation coefficients using either "
                "'pearson' (standard correlation coefficient), "
                "'correlation' (shrinkage correlation estimate), "
                "'precision' (shrinkage inverse correlation estimate), or "
                "'spearman' (Spearman rank correlation)."))
    parser.add_argument('--trials', type=int, required=True,
            help=('Number of times to run SBM from which to chose the lowest description lenght.'))
    parser.add_argument('--sigma', type=float, required=True,
            help=('The correlation cut-off to use.'))
    parser.add_argument('--data', required=True,
            help=('Path to the full graph generated by graphtool.'))
    parser.add_argument('--tissue', required=True,
            choices=('head', 'body'),
            help='Tissue being analysed.')
        
    args = parser.parse_args()

    g = load_graph(args.data)

    # Prune by treashold correlation
    tv = filterByEdge(g, args.correlation, args.sigma, True)
    g = Graph(tv, prune = True)

    N = len(g.get_vertices())
    print(N)
    Et = (N * N - N)/2
    E = len(g.get_edges())
    print(E/Et)

    plot_file = "../data/output/SBM/plots/graph_plot_" + args.tissue + "_cutoff-" + args.correlation + "_val-" + str(args.sigma) + "_non-clustered.png"
    graph_draw(g, vertex_text=g.vertex_index, edge_pen_width=g.edge_properties[args.correlation], output=plot_file)

    corr = g.edge_properties[args.correlation]

    run_non_nested_SBM(g, corr, args)
    run_nested_SBM(g, corr, args)

    