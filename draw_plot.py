from tools.tools import draw_plot, draw_DEG
import random
import numpy as np

# fix random seed
random.seed(2024)
np.random.seed(2024)

resolutions = [0.5,0.3]
datasets = ["GSE206785","GSE206785_tumor","kidney","pancreas"]
genelists = ["mt_genes", "all_genes"]
cluster_methods = ["louvain"]

draw_plot('/home/foundation/program/Foundation/record/figures/sankey-compare',datasets,genelists,cluster_methods,resolutions)


resolutions = [0.5,0.3]
datasets = ["GSE206785","GSE206785_tumor","kidney","pancreas"]
genelists = ["mt_genes(1162)"]
embedding_methods = [ "scVI"]
cluster_methods = ["louvain"]

# Call the function to process all combinations
draw_DEG("/home/foundation/program/Foundation/record/figures/DEG", datasets, genelists, cluster_methods, resolutions, embedding_methods)

