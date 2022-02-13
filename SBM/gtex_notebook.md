# gTex notebook 

## Graphs

Generated graphs are in data/output/SBM/gtex/graph, with only spearman correlations in the edge properties.

FDR trimming of the edges did not result is a substantial reduction in the graph density, unlike the fly data. Using and FDR of 1% on the SKIN graph lead to a density of 0.75, which is what I settled on.

Running the initial nested weighted SBM took 260min (about 4h20min). Now running annealing over-night