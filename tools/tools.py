import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore", category=Warning)

import scanpy as sc
import copy
from joblib import Parallel, delayed

sc.settings.set_figure_params(dpi=300, facecolor='white')

# need to make 'celltype' obs exist in adata
def print_adata_and_num_classes(data):
    """
    Print information about the AnnData object and the number of classes for each cell type.

    Parameters:
    data (AnnData): The AnnData object.

    Returns:
    None
    """
    print("adata shape:", data.shape)
    print("adata obs keys:", data.obs_keys())
    print("adata var keys:", data.var_keys())
    print("adata.X shape:", data.X.shape)
    print("adata.obs shape:", data.obs.shape)
    print("adata.var shape:", data.var.shape)
    print("Celltype:", data.obs["celltype"].unique())

    celltype_counts = data.obs["celltype"].value_counts()
    for celltype, count in celltype_counts.items():
        print(f"Celltype: {celltype}, Count: {count}")
    print(len(celltype_counts))


def display_umap(
    adata,
    color_column="celltype",
    save_path=None,
):
    """
    Display a UMAP plot for the given AnnData object.

    Parameters:
    adata (AnnData): The AnnData object containing the data.
    color_column (str): The column name in adata.obs to use for coloring the points.
    save_path (str): Path to save the plot. If None, the plot will not be saved.
    """
    p_adata = adata.copy()
    if "X_pca" not in p_adata.obsm.keys():
        print("Running PCA...")
        sc.tl.pca(p_adata, svd_solver="arpack")

    if "neighbors" not in p_adata.uns:
        print("Computing neighbors...")
        sc.pp.neighbors(p_adata, n_neighbors=10, n_pcs=40)

    if "X_umap" not in p_adata.obsm.keys():
        print("Running UMAP...")
        sc.tl.umap(p_adata)

    print("Plotting UMAP...")
    sc.pl.umap(p_adata, color=color_column, save=save_path if save_path else None)

def visualize_umap(adata, embedding_key, save_path):
    """
    Visualize an AnnData object using UMAP and save the result.
    
    Parameters:
    adata: AnnData object
    embedding_key: Key (column name) used for UMAP
    save_path: Path to save the image
    """

    # Compute neighbors and UMAP embedding
    sc.pp.neighbors(adata, use_rep=embedding_key)
    sc.tl.umap(adata)
    
    # Plot UMAP and save the image
    sc.pl.umap(adata, color="celltype", save=save_path if save_path else None, show=False)

def load_gene_list(file_path):
    """
    Load a gene list from a file.

    Parameters:
    - file_path: path to the file containing the gene list
    """
    gene_data = pd.read_json(file_path)
    if isinstance(gene_data, pd.DataFrame):
        return gene_data.columns.tolist()
    else:
        return list(gene_data.keys())
    
import pandas as pd
from pysankey2 import Sankey

def plot_sankey(data, cluster_column, celltype_column, save_path):
    """
    Plot a Sankey diagram to visualize the relationship between clusters and cell types.

    Parameters:
    - data: AnnData object
    - cluster_column: column name in data.obs containing cluster information
    - celltype_column: column name in data.obs containing cell type information
    - save_path: path to save the plot
    """
    res = pd.DataFrame({cluster_column: data.obs[cluster_column], 'celltype': data.obs[celltype_column]})
    sky = Sankey(res, colorMode="global", stripColor='left')

    fig, ax = sky.plot()
    fig.savefig(save_path)
    print(f"Sankey diagram saved to {save_path}")

from scipy.sparse import csr_matrix

def convert_chunk_to_dense(chunk):
    """
    Convert a sparse matrix chunk to a dense matrix and ensure the data type is float64.
    
    Parameters
    ----------
    chunk: csr_matrix
        The sparse matrix chunk to convert.
    
    Returns
    -------
    np.ndarray
        The dense matrix with data type as float64.
    """
    dense_chunk = chunk.toarray()
    return dense_chunk.astype(np.float64)

def convert_X_to_dense(adata):
    """
    Convert the X attribute of the AnnData object to a dense matrix and ensure the data type is float64.
    
    Parameters
    ----------
    adata: AnnData
        The AnnData object containing the X attribute.
    
    Returns
    -------
    adata: AnnData
        The modified AnnData object with X attribute as a dense matrix and data type as float64.
    """
    adata_new = adata.copy()
    if isinstance(adata.X, csr_matrix):
        n_jobs = -1  # Use all available CPU cores
        n_chunks = 5  # Number of chunks to split the matrix into
        chunk_size = adata.X.shape[0] // n_chunks
        last_chunk_size = adata.X.shape[0] - chunk_size * (n_chunks - 1)
        
        print("Splitting the matrix into multiple chunks")
        chunks = [adata.X[i*chunk_size:(i+1)*chunk_size] if i < n_chunks - 1 else adata.X[i*chunk_size:] for i in range(n_chunks)]
        
        print("Converting chunks to dense matrices in parallel")
        dense_chunks = Parallel(n_jobs=n_jobs)(delayed(convert_chunk_to_dense)(chunk) for chunk in chunks)
        
        print("Pre-allocating a large array and filling it")
        adata_new.X = np.empty((adata.X.shape[0], adata.X.shape[1]), dtype=np.float64)
        start_row = 0
        for dense_chunk in dense_chunks:
            end_row = start_row + dense_chunk.shape[0]
            adata_new.X[start_row:end_row, :] = dense_chunk
            start_row = end_row
        
        print("Conversion completed")
    else:
        print("X is already dense, if needed, convert to float64")
        adata_new.X = adata_new.X.astype(np.float64)
    
    # Print detailed information about the X attribute to confirm the changes
    print(f"Data type (dtype): {adata_new.X.dtype}")

    return adata_new

def BasicProprecess(adata,qc_min_genes=200):
    """
    do basic filtering
    :param adata: adata object
    :param qc_min_genes: filter cell, which expresses genes number less than this paramters
    :return: adata object, and plot stored in figures sub-folder
    """

    sc.pp.filter_cells(adata, min_genes=qc_min_genes)

    return adata

import anndata

def main_gene_selection(adata, gene_list):
    """
    Describe:
        Rebuild the input adata to select target genes encoding proteins.
    Parameters:
        adata->`~anndata.AnnData` object: adata with var index_name by gene symbol.
        gene_list->list: Wanted target genes.
    Returns:
        adata_new->`~anndata.AnnData` object: Modified AnnData object with selected genes.
        to_fill_columns->list: Zero padding genes.
    """

    # Extract the expression matrix and the genes
    X_df = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)
    
    # Identify genes that need to be added (zero-padded)
    to_fill_columns = list(set(gene_list) - set(X_df.columns))
    
    # Create a DataFrame with zero padding for missing genes
    padding_df = pd.DataFrame(np.zeros((X_df.shape[0], len(to_fill_columns))), 
                              columns=to_fill_columns, 
                              index=X_df.index)
    
    # Concatenate the original DataFrame with the zero-padded DataFrame
    X_df = pd.concat([X_df, padding_df], axis=1)
    
    # Select the genes in the order specified by gene_list
    X_df = X_df[gene_list]
    
    # Create a new AnnData object with the modified expression matrix
    adata_new = anndata.AnnData(X=X_df.values, obs=adata.obs, var=pd.DataFrame(index=gene_list))
    
    return adata_new, to_fill_columns


def plot_sankey2(data1, data2, cluster_column, save_path):
    res = pd.DataFrame({
        f'{cluster_column}_1': data1.obs[cluster_column],
        f'{cluster_column}_2': data2.obs[cluster_column]
    })
    res.columns = ['source', 'target']

    sky = Sankey(res, colorMode="global", stripColor='left')

    fig, ax = sky.plot()
    fig.savefig(save_path)
    print(f"Sankey diagram saved to {save_path}")


import os
def draw_plot(save_path,datasets,genelists,cluster_methods,resolutions):
    h5ad_directory = "/home/foundation/program/foundation-new/record/temp-h5ad"
    sc.settings.figdir = "/home/foundation/program/foundation-new/record/figures/UMAP"
    sc.settings.plot_prefix = ""
    for cluster_method in cluster_methods:
        for dataset in datasets:
            for genelist in genelists:
                for resolution in resolutions:
                    # scVI -------------------------------
                    file_name = f"{dataset}-{genelist}-{cluster_method}-{resolution}-scVI.h5ad"
                    file_path = os.path.join(h5ad_directory, file_name)

                    if not os.path.exists(file_path):
                        ValueError(f"File not found: {file_path}")

                    model_data = sc.read_h5ad(file_path)
                    
                    sc.pp.neighbors(model_data, use_rep="X_scVI")
                    sc.tl.umap(model_data)

                    # draw two UMAP plots and save (one for cluster and one for cell type)
                    sc.pl.umap(model_data, color=[cluster_method], title=None, save=f"{dataset}-{genelist}-{cluster_method}-{resolution}-scVI-Cluster.pdf")
                    sc.pl.umap(model_data, color=['celltype'], title=None, save=f"{dataset}-{genelist}-{cluster_method}-{resolution}-scVI-Celltype.pdf")
                    
                    # draw sankey diagram and save
                    plot_sankey(model_data, cluster_column= cluster_method, celltype_column='celltype', save_path=f"/home/foundation/program/foundation-new/record/figures/sankey/{dataset}-{genelist}-{cluster_method}-{resolution}-scVI.pdf")


                    # PCA --------------------------------
                    file_name = f"{dataset}-{genelist}-{cluster_method}-{resolution}-PCA.h5ad"
                    file_path = os.path.join(h5ad_directory, file_name)

                    if not os.path.exists(file_path):
                        ValueError(f"File not found: {file_path}")

                    PCA_data = sc.read_h5ad(file_path)

                    sc.pp.neighbors(PCA_data, use_rep="X_pca")
                    sc.tl.umap(PCA_data)

                    # draw two UMAP plots and save (one for cluster and one for cell type)
                    sc.pl.umap(PCA_data, color=[cluster_method], title=None, save=f"{dataset}-{genelist}-{cluster_method}-{resolution}-PCA-Cluster.pdf")
                    sc.pl.umap(PCA_data, color=['celltype'], title=None, save=f"{dataset}-{genelist}-{cluster_method}-{resolution}-PCA-Celltype.pdf")
                    
                    # draw sankey diagram and save
                    plot_sankey(PCA_data, cluster_column= cluster_method, celltype_column='celltype', save_path=f"/home/foundation/program/foundation-new/record/figures/sankey/{dataset}-{genelist}-{cluster_method}-{resolution}-PCA.pdf")


                    # scGPT --------------------------------
                    file_name = f"{dataset}-{genelist}-{cluster_method}-{resolution}-scGPT.h5ad"
                    file_path = os.path.join(h5ad_directory, file_name)

                    if not os.path.exists(file_path):
                        ValueError(f"File not found: {file_path}")

                    scGPT = sc.read_h5ad(file_path)

                    sc.pp.neighbors(scGPT, use_rep="X_scGPT")
                    sc.tl.umap(scGPT)

                    # draw two UMAP plots and save (one for cluster and one for cell type)
                    sc.pl.umap(scGPT, color=[cluster_method], title=None, save=f"{dataset}-{genelist}-{cluster_method}-{resolution}-scGPT-Cluster.pdf")
                    sc.pl.umap(scGPT, color=['celltype'], title=None, save=f"{dataset}-{genelist}-{cluster_method}-{resolution}-scGPT-Celltype.pdf")
                    
                    # draw sankey diagram and save
                    plot_sankey(scGPT, cluster_column= cluster_method, celltype_column='celltype', save_path=f"/home/foundation/program/foundation-new/record/figures/sankey/{dataset}-{genelist}-{cluster_method}-{resolution}-scGPT.pdf")


                    # scGPT_allgene --------------------------------
                    file_name = f"{dataset}-{genelist}-{cluster_method}-{resolution}-scGPT_allgene.h5ad"
                    file_path = os.path.join(h5ad_directory, file_name)

                    if not os.path.exists(file_path):
                        ValueError(f"File not found: {file_path}")

                    scGPTallgene = sc.read_h5ad(file_path)

                    sc.pp.neighbors(scGPTallgene, use_rep="X_scGPT")
                    sc.tl.umap(scGPTallgene)

                    # draw two UMAP plots and save (one for cluster and one for cell type)
                    sc.pl.umap(scGPTallgene, color=[cluster_method], title=None, save=f"{dataset}-{genelist}-{cluster_method}-{resolution}-scGPT_B-Cluster.pdf")
                    sc.pl.umap(scGPTallgene, color=['celltype'], title=None, save=f"{dataset}-{genelist}-{cluster_method}-{resolution}-scGPT_B-Celltype.pdf")
                    
                    # draw sankey diagram and save
                    plot_sankey(scGPTallgene, cluster_column= cluster_method, celltype_column='celltype', save_path=f"/home/foundation/program/foundation-new/record/figures/sankey/{dataset}-{genelist}-{cluster_method}-{resolution}-scGPT_B.pdf")


                    # sankey (model VS PCA, model VS scGPT, model VS scGPT_allgene)
                    plot_sankey2(model_data, PCA_data, cluster_column=cluster_method, save_path=os.path.join(save_path, f"modelVSPCA-{dataset}-{genelist}-{cluster_method}-{resolution}.pdf"))
                    plot_sankey2(model_data, scGPT, cluster_column=cluster_method, save_path=os.path.join(save_path, f"modelVSscGPT-{dataset}-{genelist}-{cluster_method}-{resolution}.pdf"))
                    plot_sankey2(model_data, scGPTallgene, cluster_column=cluster_method, save_path=os.path.join(save_path, f"modelVSscGPTallgene-{dataset}-{genelist}-{cluster_method}-{resolution}.pdf"))


# Function to rank and plot genes
def rank_and_plot_genes(adata, save_path, dataset, genelist, cluster_method, resolution, embedding_method, embedding_key, method="wilcoxon", n_genes1=20, n_genes2=5):

    # Set the save directory for scanpy plots
    sc.settings.figdir = save_path
    sc.settings.plot_prefix = ""

    # Verify if the embedding key exists in adata.obs
    if embedding_key not in adata.obs.columns:
        ValueError(f"{embedding_key} does not exist in adata.obs columns. Available keys are: {adata.obs.columns.tolist()}")
        

    # # Rank genes groups and plot the first set of genes
    # sc.tl.rank_genes_groups(adata, embedding_key, method=method)           
    # sc.pl.rank_genes_groups(adata, n_genes=n_genes1, sharey=False, save=f"{dataset}-{genelist}-{cluster_method}-{resolution}-{embedding_method}-rankgenes.pdf")

    # Rank genes groups and plot the second set of genes
    sc.tl.rank_genes_groups(adata, embedding_key, method=method)
    # sc.pl.rank_genes_groups(adata, n_genes=n_genes2, sharey=False)

    # Extract top genes
    rank_genes_groups = adata.uns["rank_genes_groups"]
    names = pd.DataFrame(rank_genes_groups['names'])
    classes = names.columns

    top_genes = {}
    for classs in classes:
        top_genes[classs] = names[classs].values[:n_genes2]

    top_genes_df = pd.DataFrame(top_genes)

    # Pick unique genes
    unique_top_genes_flat_list = []
    seen_genes = set()

    for col in top_genes_df.columns:
        for gene in top_genes_df[col]:
            if gene not in seen_genes:
                unique_top_genes_flat_list.append(gene)
                seen_genes.add(gene)

    # Plot dotplot and stacked violin plot
    sc.pl.dotplot(adata, unique_top_genes_flat_list, groupby=embedding_key, save=f"{dataset}-{genelist}-{cluster_method}-{resolution}-{embedding_method}-dotplot.pdf")
    sc.pl.stacked_violin(adata, unique_top_genes_flat_list, groupby=embedding_key, save=f"{dataset}-{genelist}-{cluster_method}-{resolution}-{embedding_method}-violinplot.pdf")

    ranked_genes_df = sc.get.rank_genes_groups_df(adata, group=None, pval_cutoff=0.05)

    top_genes_df_pvals = (ranked_genes_df
                    .sort_values(by=["group", "pvals"])
                    .groupby("group")
                    .head(15)
                    .sort_values(by=["group", "pvals"]))

    top_genes_df_logfc = (ranked_genes_df
                .sort_values(by=["group", "logfoldchanges"], ascending=[True, False])
                .groupby("group")
                .head(15)
                .sort_values(by=["group", "logfoldchanges"], ascending=[True, False]))

    return top_genes_df_pvals, top_genes_df_logfc

# Function to process all combinations of parameters
def draw_DEG(save_path, datasets, genelists, cluster_methods, resolutions, embedding_methods):
    h5ad_directory = "/home/foundation/program/foundation-new/record/temp-h5ad"
    for dataset in datasets:
        for genelist in genelists:
            for embedding_method in embedding_methods:
                for cluster_method in cluster_methods:
                    for resolution in resolutions:  
                        file_name = f"{dataset}-{genelist}-{cluster_method}-{resolution}-{embedding_method}.h5ad"
                        file_path = os.path.join(h5ad_directory, file_name)

                        if not os.path.exists(file_path):
                            ValueError(f"File not found: {file_path}")

                        adata = sc.read_h5ad(file_path)
                        top_genes_df_pvals, top_genes_df_logfc = rank_and_plot_genes(adata, save_path, dataset, genelist, cluster_method, resolution, embedding_method, embedding_key=cluster_method, method="wilcoxon", n_genes1=20, n_genes2=5)

                        top_genes_df_pvals.to_csv(os.path.join(f"/home/foundation/program/foundation-new/record/pvals/{dataset}-{genelist}-{cluster_method}-{resolution}-{embedding_method}-topgenes-pvals.csv"))
                        top_genes_df_logfc.to_csv(os.path.join(f"/home/foundation/program/foundation-new/record/logfc/{dataset}-{genelist}-{cluster_method}-{resolution}-{embedding_method}-topgenes-logfc.csv"))