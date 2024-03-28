import os
os.environ["OMP_NUM_THREADS"] = "4"

from graph_tool.all import *
import pandas as pd
import numpy as np
import matplotlib.cm as mpl

from trim_networks import *
from fit_sbm import *

def FilterByGraph(g, g_ref):
    genes = g.vp.genes

    keep = g.new_vertex_property("bool", True)
    for v in g.get_vertices():
        if genes[v] not in g_ref.vp.genes:
            keep[v] = False
    tv = GraphView(g, vfilt=keep)
    g = Graph(tv, prune = True)
    return g


def KeepRandomVertices(g, num_vertices_to_keep):
    keep = g.new_vertex_property("bool")
    vertices = list(g.get_vertices())
    
    # Select a random subset of vertices
    subset_vertices = np.random.choice(vertices, num_vertices_to_keep, replace=False) 
    
    # Mark the selected vertices to keep
    for v in subset_vertices:
        keep[v] = True
    
    # Create a new graph view with the filtered vertices
    tv = GraphView(g, vfilt=keep)
    
    # Create a new graph from the graph view and prune unused vertices
    g_subset = Graph(tv, prune=True)
    
    return g_subset


def KeepRandomSubset(g, subset_ratio):
    keep = g.new_vertex_property("bool")
    vertices = list(g.get_vertices())
    num_vertices = len(vertices)
    num_subset = int(num_vertices * subset_ratio)
    
    # Select a random subset of vertices
    subset_vertices = np.random.choice(vertices, num_subset, replace=False)
    
    # Mark the selected vertices to keep
    for v in subset_vertices:
        keep[v] = True
    
    # Create a new graph view with the filtered vertices
    tv = GraphView(g, vfilt=keep)
    
    # Create a new graph from the graph view and prune unused vertices
    g_subset = Graph(tv, prune=True)
    
    return g_subset

def density(g):
    N = len(g.get_vertices())
    Et = (N * N - N)/2
    E = len(g.get_edges())
    density = E/Et
    return density

folder = "../data/output/SBM/graphs/"
fdr_list = [0.5, 0.25, 1e-1, 5e-2, 1e-2]
#fdr_list = [0.25, 1e-1]
iter = 1000
nodes = 300
seed = 0

#g_body = load_graph(folder + "VOOMCounts_CPM1_body_ctrl_249ind_counts3M_covfree_Aug3121.xml.gz")

output_folder = "FDRplots/" + "nodes-" + str(nodes) + "/" + "seed-" + str(seed) + "/"
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
plot_folder = os.path.join(output_folder, "plots/")
if not os.path.exists(plot_folder):
    os.makedirs(plot_folder)
graph_folder = os.path.join(output_folder, "graphs/")
if not os.path.exists(graph_folder):
    os.makedirs(graph_folder)

# Check if the file exists
if os.path.exists(output_folder + "SBM_head_" + str(nodes) + "_vertices.xml.gz"):
    # If it does, load the graph from the file
    print("Loading small head graph")
    g_total = load_graph(output_folder + "SBM_head_" + str(nodes) + "_vertices.xml.gz")
else:
    print("Loading head graph")
    g_head = load_graph(folder + "VOOMCounts_CPM1_head_ctrl_248ind_counts3M_covfree_Aug3121.xml.gz")
    # Filter the graph to keep only the top 150 vertices
    print("Filtering graph vertices")
    np.random.seed(seed)
    g_total = KeepRandomVertices(g_head, nodes)
    g_total.save(output_folder + "SBM_head_" + str(nodes) + "_vertices.xml.gz")

corr = g_total.edge_properties["spearman"]
g_total.ep.weight = g_total.new_edge_property("double", (2*np.arctanh(corr.a)))

print("Filtering graph by FDR")
g = {}
for fdr in fdr_list:
    g[fdr] = filterByFDR(g_total, fdr, True)

print("Filtering by smallest graph")
for i in fdr_list:
    g[i] = FilterByGraph(g[i], g[fdr_list[-1]])

print("Calculating layout")
pos = arf_layout(g[fdr_list[-1]])

state = {}
for i in fdr_list:

    if os.path.exists(graph_folder + "SBM_FDR-" + str(i) +  "_blocks.dill"):
        print("Loading block state for FDR: ", i)
        with open(graph_folder + "SBM_FDR-" + str(i) +  "_blocks.dill", 'rb') as fh:
            block_state = dill.load(fh)
        state[i] = minimize_nested_blockmodel_dl(g[i], init_bs = block_state,
                                                 state_args=dict(recs=[g[i].ep.weight],
                                                 rec_types=["real-normal"]))
    else:
        print("Fitting SBM for FDR: ", i)
        state[i] = minimize_nested_blockmodel_dl(g[i], 
                                                state_args=dict(recs=[g[i].ep.weight],
                                                rec_types=["real-normal"]))
        
        print("Annealing for FDR: ", i)
        mcmc_anneal(state[i], beta_range=(1, 10), niter=iter,
                    mcmc_equilibrate_args=dict(force_niter=100))

        print("Equilibrating for FDR: ", i)
        mcmc_equilibrate(state[i], wait=iter, nbreaks=1, mcmc_args=dict(niter=10))

    # collect nested partitions
    bs = []

    def collect_partitions(s):
        global bs
        bs.append(s.get_bs())

    print("Collecting MCMC marginals for FDR: ", i)
    mcmc_equilibrate(state[i], force_niter=iter*10, mcmc_args=dict(niter=10),
                     callback=collect_partitions)

    print("Disambiguating partitions for FDR: ", i)
    pmode = PartitionModeState(bs, nested=True, converge=True)
    pv = pmode.get_marginal(g[i])

    g[i].vp.pv = g[i].new_vertex_property("int")
    g[i].vp.pv = pv

    print("Getting consensus estimate for FDR: ", i)
    bs = pmode.get_max_nested()
    state[i] = state[i].copy(bs=bs)

    print("Drawing graph for FDR: ", i)
    g[i].vp.level0 = g[i].new_vertex_property("int", np.array(state[i].get_bs()[0]))
    state[i].draw(output=plot_folder + "circle-SBM_FDR-" + str(i) + "_graph.png", 
                  eorder=g[i].ep.weight,
                  edge_color=prop_to_size(g[i].ep.weight, mi=-4, ma=4, power=1, log=False),
                  ecmap=(mpl.inferno, .6), 
                  edge_gradient=[], 
                  hvertex_size = 15,
                  hedge_pen_width = 3,
                  vertex_color = g[i].vp.level0,
                  vertex_fill_color = g[i].vp.level0)

    print("Drawing graph for FDR: ", i)
    graph_draw(g[i], 
               pos = pos, 
               edge_pen_width = g[i].ep.weight, 
               eorder=g[i].ep.weight,
               vertex_shape="pie", vertex_pie_fractions=pv,
               edge_color=prop_to_size(g[i].ep.weight, mi=-4, ma=4, power=1, log=False),
               ecmap=(mpl.inferno, .6), 
               edge_gradient=[], 
               output = plot_folder + "network-SBM_FDR-" + str(i) + "_graph.png")

    print("Saving block state for FDR: ", i)
    block_state = state[i].get_bs()
    output_file = graph_folder + "SBM_FDR-" + str(i) +  "_blocks.dill"
    with open(output_file, 'wb') as fh:
        dill.dump(block_state, fh, recurse=True)

    print("Saving graph for FDR: ", i)
    g[i].save(graph_folder + "SBM_FDR-" + str(i) +  "_graph.xml.gz")

n_block = {}
for i in fdr_list:
    bs = state[i].get_bs()[0]
    # get number of distinct blocks in bs
    n_block[i] = len(set(bs))
# get key of shortest block
key = min(n_block, key=n_block.get)

print("Aligning partitions")
ref = state[key].get_bs()[0]
aligned_bs = {} 
for i in fdr_list:
    bs = state[i].get_bs()[0]
    aligned_bs[i] = align_partition_labels(bs, ref)

for i in fdr_list:
    print("Drawing graph for FDR: ", i)
    g[i].vp.aligned_bs = g[i].new_vertex_property("int", aligned_bs[i])
    graph_draw(g[i], 
               pos = pos, 
               edge_pen_width = g[i].ep.weight, 
               eorder=g[i].ep.weight,
               vertex_color=g[i].vp.aligned_bs, vertex_fill_color=g[i].vp.aligned_bs,
               edge_color=prop_to_size(g[i].ep.weight, mi=-4, ma=4, power=1, log=False),
               ecmap=(mpl.inferno, .6), 
               edge_gradient=[], 
               output = plot_folder + "network-SBM_FDR-" + str(i) + "_graph.png")
    
    graph_draw(g[i], 
               pos = pos, 
               vertex_shape="pie", vertex_pie_fractions=g[i].vp.pv,
               edge_pen_width = g[i].ep.weight, 
               eorder=g[i].ep.weight,
               edge_color=prop_to_size(g[i].ep.weight, mi=-4, ma=4, power=1, log=False),
               ecmap=(mpl.inferno, .6), 
               edge_gradient=[], 
               output = plot_folder + "network-SBM_FDR-" + str(i) + "_graph_pie.png")
    
    print("Drawing circle graphs for FDR: ", i)
    state[i].draw(output=plot_folder + "circle-SBM_FDR-" + str(i) + "_graph.png", 
                  eorder=g[i].ep.weight,
                  edge_color=prop_to_size(g[i].ep.weight, mi=-4, ma=4, power=1, log=False),
                  ecmap=(mpl.inferno, .6), 
                  edge_gradient=[], 
                  hvertex_size = 15,
                  hedge_pen_width = 3,
                  vertex_color = g[i].vp.aligned_bs,
                  vertex_fill_color = g[i].vp.aligned_bs)

print("Done")