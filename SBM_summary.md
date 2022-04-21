# DC-SBM (degree corrected stochastic block model)

Finds the network partition that best describes the network (shortest description lengh, less bits used to describe the same network, the partition that conveis the largest amount of information about the network). This is done by grouping nodes according to the probability that we observe an edge between and within blocks. There is no requiremennt that blocks be more cohessive or assortative, only that the nodes in the blocks have a similar pattern of connections with other blocks (and with themselves, but less so in the degree connected model). In the degree-corrected model there are node-level parameters that account for the degree distribution inside each block, so nodes that have small and large degrees can be clustered together in the same block. The advantage is that this will pickup on whatever pattern in the network that is usefull in describing it (but it's up to us to interpret the blocks and see if they are interesting and in what way). We could end up with a mixture of different types of block: assortative modules, bridges between modules, other things I can't think of... 

# MMC 

Modularity maximization using a transformed similarity matrix to increase the constrast between high similarity and low similarity. I suspect this could be modified to work reasonably well with the Assortative-SBM with little work. But it's just modularity maximization (with all the usual under- and over- fitting problems)

# Assortative SBM (also called Planted Partition model)

Similar to the block model, but with the additional constraint that the number of edges inside a block is expected to be larger than the number os edges between blocks. This looks for traditional "modular" partition, and is an improved version of the usual modularity maximization (but without over-fitting (having too many modules) and under-fitting (having too few modules due to the resolution limit)). This is more restrictive than the DC-SBM.

# WGCNA

Basically hierarchical clustering and tresholding the resulting tree at some arbritrary height. 

In the expression data it ends up a lot of singletons or unconnected nodes because there is a large number of genes with low correlation with everything else that are not close enough to other groups to be bundled together. I predict all that "fuzz" around the main network is just getting discarted (but could be easily placed in a module by the other methods since it's just connected to a single cohesive module).

