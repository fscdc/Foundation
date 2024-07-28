from tools.tools import visualize_umap, plot_sankey, analyze_datasets
import random
import numpy as np

# fix random seed
random.seed(2024)
np.random.seed(2024)

resolutions = [0.5,0.3]
datasets = ["GSE206785","GSE206785_tumor","kidney","pancreas"]
genelists = ["mt_genes", "all_genes"]
embedding_methods = ["scVI"]
cluster_methods = ["louvain"]

# 调用函数
analyze_datasets(datasets,genelists,cluster_methods,resolutions,embedding_methods)


import scanpy as sc
import pandas as pd
from pysankey2 import Sankey
import os



def plot_sankey(data1, data2, cluster_column, save_path):
    res = pd.DataFrame({
        f'{cluster_column}_1': data1.obs[cluster_column],
        f'{cluster_column}_2': data2.obs[cluster_column]
    })
    res.columns = ['source', 'target']

    sky = Sankey(res, colorMode="global", stripColor='left')

    fig, ax = sky.plot()
    fig.savefig(save_path)
    print(f"Sankey diagram saved to {save_path}")

def process_h5ad_files(h5ad_directory, save_path,datasets,genelists,cluster_methods,resolutions):
    for cluster_method in cluster_methods:
        for dataset in datasets:
            for genelist in genelists:
                for resolution in resolutions:
                    file_name = f"{dataset}-{genelist}-{cluster_method}-{resolution}-scVI.h5ad"
                    file_path = os.path.join(h5ad_directory, file_name)
                    model = sc.read_h5ad(file_path)
                    
                    file_name = f"{dataset}-{genelist}-{cluster_method}-{resolution}-PCA.h5ad"
                    file_path = os.path.join(h5ad_directory, file_name)
                    PCA = sc.read_h5ad(file_path)

                    file_name = f"{dataset}-{genelist}-{cluster_method}-{resolution}-scGPT.h5ad"
                    file_path = os.path.join(h5ad_directory, file_name)
                    scGPT = sc.read_h5ad(file_path)

                    file_name = f"{dataset}-{genelist}-{cluster_method}-{resolution}-scGPT_allgene.h5ad"
                    file_path = os.path.join(h5ad_directory, file_name)
                    scGPTallgene = sc.read_h5ad(file_path)

                    # Call the plot_sankey function
                    plot_sankey(model, PCA, cluster_column=cluster_method, save_path=os.path.join(save_path, f"modelVSPCA-{dataset}-{genelist}-{cluster_method}-{resolution}.pdf"))
                    plot_sankey(model, scGPT, cluster_column=cluster_method, save_path=os.path.join(save_path, f"modelVSscGPT-{dataset}-{genelist}-{cluster_method}-{resolution}.pdf"))
                    plot_sankey(model, scGPTallgene, cluster_column=cluster_method, save_path=os.path.join(save_path, f"modelVSscGPTallgene-{dataset}-{genelist}-{cluster_method}-{resolution}.pdf"))

# Input parameters
# resolutions = [0.3, 0.5]
# datasets = ["GSE206785", "GSE261157", "GSE206785_tumor"]
# genelists = ["all_genes", "mt_genes_1(1683)", "mt_genes_2(1391)"]
# cluster_methods = ["louvain", "leiden"]


resolutions = [0.3,0.5]
datasets = ["pancreas"]
genelists = ["mt_genes_2(1391)","mt_genes_1(1683)"]
embedding_methods = ["PCA", "scGPT", "scGPT_allgene", "scVI"]
cluster_methods = ["louvain"]

process_h5ad_files('/home/foundation/program/foundation/temp-h5ad','/home/foundation/program/foundation/plot/sankey-compare',datasets,genelists,cluster_methods,resolutions)


# visualize_umap(
#     adata_new, 
#     embedding_key="X_scVI", 
#     save_path=f"{args.dataset}-{args.gene_list}-{args.cluster_method}-{args.resolution}-scVI.pdf"
# )