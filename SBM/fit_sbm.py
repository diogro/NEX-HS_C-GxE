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

def create_Block_df(g, corr, state):
    genes = g.vertex_properties["genes"]
    block_df = pd.DataFrame(columns=('Gene', "Degree", "E_corr", 'Block'))
    for v in g.vertex_index:
        line = [genes[v]]
        line.append(g.get_total_degrees([v])[0])
        line.append(np.mean(np.abs(g.get_all_edges(v, [corr] )[:,2])))
        line.append(state.get_blocks()[v])
        block_df.loc[v] = line
    return block_df

def get_group(x, state):
    levels = state.get_levels()
    n_levels = 7#len(levels)
    r = np.zeros(n_levels)
    r[0] = levels[0].get_blocks()[x]
    for i in range(1, n_levels):
        r[i] = levels[i].get_blocks()[r[i-1]]
    r = r.astype(int)
    return r
def create_nestedBlock_df(g, corr, state):
    genes = g.vertex_properties["genes"]
    nested_block_df = pd.DataFrame(columns=('Gene', "Degree", "E_corr", 'B1', "B2", "B3", "B4", "B5", "B6", "B7"))
    for v in g.vertex_index:
        line = [genes[v]]
        line.append(g.get_total_degrees([v])[0])
        line.append(np.mean(np.abs(g.get_all_edges(v, [corr] )[:,2])))
        [line.append(i) for i in get_group(v, state)]
        nested_block_df.loc[v] = line
    return nested_block_df

def run_non_nested_SBM(g, corr, args):
    print("Running non-nested model.")
    state_min_nn = minimize_blockmodel_dl(g)
    mcmc_equilibrate(state_min_nn, wait=args.wait, mcmc_args=dict(niter=10))

    plot_file = out_folder + "plots/graph_plot_" + out_label + "_clustered-non-hierarchical-SBM.png"
    state_min_nn.draw(output = plot_file)

    block_df = create_Block_df(g, corr, state_min_nn)
    output_file = out_folder + "clustering/" + out_label + "_non-hierarchical-SBM.csv"
    block_df.to_csv(output_file)

def run_nested_SBM(g, corr, args, blocks=None):

    print("Creating nested model...")
    if args.layer is True:
        state_min = minimize_nested_blockmodel_dl(g, init_bs=blocks, 
                                                  state_args=dict(base_type=LayeredBlockState,
                                                                  ec=g.ep.layer,
                                                                  layers=True,
                                                                  recs=[g.ep.z_s],
                                                                  rec_types=["real-normal"]))
    else:
        state_min = minimize_nested_blockmodel_dl(g, init_bs=blocks, 
                                                  state_args=dict(recs=[g.ep.z_s],
                                                                  rec_types=["real-normal"]))


    initial_entropy = state_min.entropy()
    print("Created nested model. Initial entropy: " + str(initial_entropy))

    if args.annealing > 0:
        S1 = state_min.entropy()
        print("Starting annealing...")
        mcmc_anneal(state_min, beta_range=(1, 10), niter=args.annealing, 
                    mcmc_equilibrate_args=dict(force_niter=10))
        S2 = state_min.entropy()
        print("Improvement from annealing:", S2 - S1)
        print("Final entropy after annealing: " + str(state_min.entropy()))

        if state_min.entropy() < initial_entropy:
            plot_file = out_folder + "plots/graph_plot_" + out_label + "_clustered-hierarchical-SBM.png"
            state_min.draw(output = plot_file)

            nested_block_df = create_nestedBlock_df(g, corr, state_min)
            output_file = out_folder + "clustering/" + out_label + "_hierarchical-SBM.csv"
            nested_block_df.to_csv(output_file)

            block_state = state_min.get_bs()
            output_file = out_folder + "clustering/" + out_label + "_hierarchical-SBM.dill"
            with open(output_file, 'wb') as fh:
                dill.dump(block_state, fh, recurse=True)


    if args.wait > 0:
        print("Starting MCMC equilibration...")
        initial_entropy = state_min.entropy()
        mcmc_equilibrate(state_min, wait=args.wait, mcmc_args=dict(niter=10))
        print("Improvement from equilibration:", state_min.entropy() - initial_entropy)
        print("Final entropy after equilibration: " + str(state_min.entropy()))

        if state_min.entropy() < initial_entropy:
            plot_file = out_folder + "plots/graph_plot_" + out_label + "_clustered-hierarchical-SBM.png"
            state_min.draw(output = plot_file)

            nested_block_df = create_nestedBlock_df(g, corr, state_min)
            output_file = out_folder + "clustering/" + out_label + "_hierarchical-SBM.csv"
            nested_block_df.to_csv(output_file)

            block_state = state_min.get_bs()
            output_file = out_folder + "clustering/" + out_label + "_hierarchical-SBM.dill"
            with open(output_file, 'wb') as fh:
                dill.dump(block_state, fh, recurse=True)

    if args.mcmc == True:
        bs = []
        h = [np.zeros(g.num_vertices() + 1) for s in state_min.get_levels()]

        def collect_partitions(s):
            bs.append(s.get_bs())
            for l, sl in enumerate(s.get_levels()):
                B = sl.get_nonempty_B()
                h[l][B] += 1

        print("Starting MCMC.")
        # Now we collect the marginals for exactly args.iter*10 sweeps
        mcmc_equilibrate(state_min, force_niter=args.iter, mcmc_args=dict(niter=10),
                         callback=collect_partitions)

        pmode = PartitionModeState(bs, nested=True, converge=True)
        pv = pmode.get_marginal(g)

        # Get consensus estimate
        bs_max = pmode.get_max_nested()

        state = state_min.copy(bs=bs_max)

        print("Description lenght improvement in mcmc: " + str(state_min.entropy() - state.entropy()))
        print("Final entropy after MCMC: " + str(state.entropy()))

        plot_file = out_folder + "plots/graph_plot_" + out_label + "_mcmc_mode_clustered-hierarchical-SBM.png"
        state.draw(output = plot_file)

        nested_block_df = create_nestedBlock_df(g, corr, state)
        output_file = out_folder + "clustering/" + out_label + "_mcmc_mode_hierarchical-SBM.csv"
        nested_block_df.to_csv(output_file)

        for i in range(7):
            B = state.get_levels()[i].get_nonempty_B()
            e_mat = state.get_levels()[i].get_matrix().todense()
            output_file = out_folder + "clustering/" + out_label + "_mcmc_mode_hierarchical-SBM_e_matrix_level" + str(i) + ".csv"
            pd.DataFrame(e_mat).to_csv(output_file)

        output_file = out_folder + "clustering/" + out_label + "_mcmc_mode_hierarchical-SBM_levels_histogram.csv"
        pd.DataFrame(h).to_csv(output_file)

        block_state = state.get_bs()
        output_file = out_folder + "clustering/" + out_label + "_mcmc_mode_hierarchical-SBM.dill"
        with open(output_file, 'wb') as fh:
            dill.dump(block_state, fh, recurse=True)

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
    parser.add_argument('--wait', default = 100, type=int,
            help=('Number of runs without improvement of block partition in equilibration.'))
    parser.add_argument('--graph',
            help=('Path to the full input graph generated by graph-tool.'))
    parser.add_argument('--block',
            help=('Path to the pickled infered block structure.'))
    parser.add_argument('--output', required=True,
            help=('Output label.'))
    parser.add_argument('--tissue', required=True,
            choices=('head', 'body'),
            help='Tissue being analysed.')
    parser.add_argument('--type',
            choices=('all', 'nested', 'non-nested', 'none'), default="nested",
            help='Type of SBM model.  Default nested.')
    parser.add_argument('--layer', dest='layer', action='store_true')
    parser.add_argument('--no-layer', dest='layer', action='store_false')
    parser.set_defaults(layer=False)
    parser.add_argument('--mcmc', dest='mcmc', action='store_true')
    parser.add_argument('--no-mcmc', dest='mcmc', action='store_false')
    parser.set_defaults(mcmc=False)
    parser.add_argument('--mcmc-iter', dest='iter', default = 10000, type=int,
    help=('Number of MCMC iterations.'))
    parser.add_argument('--annealing', dest='annealing', default = 0, type=int,
    help=('Number of MCMC simulated annealing iterations.'))
    args = parser.parse_args()

    print("Loading graph...")
    g = load_graph(args.graph)
    if args.block is None:
        bs = None
    else:
        print("Loading blocks...")
        with open (args.block, "rb") as fh:
            bs = dill.load(fh)

    corr = g.edge_properties[args.correlation]
    g.ep.positive = g.new_edge_property("int", (np.sign(corr.a) + 1)/2)
    g.ep.layer = g.new_edge_property("int16_t", np.sign(corr.a).astype(np.int16))
    g.ep.layer.a = np.sign(corr.a).astype(np.int16)
    g.ep.z_s = g.new_edge_property("double", (2*np.arctanh(corr.a)))

    N = len(g.get_vertices())
    print(N)
    Et = (N * N - N)/2
    E = len(g.get_edges())
    print(E/Et)

    print("Min correlation: ", min(g.edge_properties[args.correlation].a))
    print("Max correlation: ", max(g.edge_properties[args.correlation].a))

    out_folder = "../data/output/SBM/"
    out_label = args.tissue + "_weights-" + args.correlation + "_" + args.output
    if args.layer is True:
        out_label = out_label + "_layered"

    print("Out file label: " + out_label)

    plot_file = out_folder + "plots/" + out_label + "_non-clustered.png"
    corr.a = np.abs(corr.a)
    if args.type != 'none':
        graph_draw(g, vertex_text=g.vertex_index, edge_pen_width=corr, output=plot_file)

    corr = g.edge_properties[args.correlation]

    if args.type == "all":
        run_non_nested_SBM(g, corr, args)
        run_nested_SBM(g, corr, args, blocks = bs)
    if args.type == 'nested':
        run_nested_SBM(g, corr, args, blocks = bs)
    if args.type == 'non-nested':
        run_non_nested_SBM(g, corr, args)
    if args.type == 'none':
        print("Test run, did nothing.")
