import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore", category=Warning)

import scanpy as sc
import copy
from joblib import Parallel, delayed

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

sc.settings.set_figure_params(dpi=300, facecolor='white')

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


def analyze_datasets(datasets, genelists,cluster_methods,resolutions,embedding_methods):
    """
    遍历数据集、基因列表、嵌入方法和分辨率，进行UMAP绘制和桑基图生成，并保存处理后的数据。
    
    Parameters
    ----------
    resolutions : list of float
        要使用的分辨率列表。
    datasets : list of str
        要处理的数据集列表。
    genelists : list of str
        要使用的基因列表。
    embedding_methods : list of str
        要使用的嵌入方法列表。
    cluster_methods : list of str
        要使用的聚类方法列表。
    """
    # 遍历数据集、基因列表、嵌入方法和分辨率

    import os
    for cluster_method in cluster_methods:
        for dataset in datasets:
            for genelist in genelists:
                for embedding_method in embedding_methods:
                    for resolution in resolutions:

                        file_name = f"{dataset}-{genelist}-{cluster_method}-{resolution}-{embedding_method}.h5ad"
                        file_path = os.path.join("/home/foundation/program/foundation-new/record/temp-h5ad/", file_name)

                        if not os.path.exists(file_path):
                            print(f"File not found: {file_path}")
                            continue

                        adata = sc.read(file_path)
                        

                        sc.settings.figdir = "/home/foundation/program/foundation-new/record/figures/umap-celltype-cluster"

                        if embedding_method == "PCA":
                            use_rep = 'X_pca'
                        elif embedding_method in {"scGPT", "scGPT_allgene"}:
                            use_rep = 'X_scGPT'
                        elif embedding_method == "scVI":
                            use_rep = 'X_scVI'

                        sc.pp.neighbors(adata, use_rep=use_rep)
                        sc.tl.umap(adata)

                        # draw two UMAP plots and save (one for cluster and one for cell type)
                        sc.pl.umap(adata, color=[cluster_method], title=f'{cluster_method} (resolution={resolution})', save=f"_cluster-{dataset}-{genelist}-{cluster_method}-{resolution}-{embedding_method}.pdf")
                        sc.pl.umap(adata, color=['celltype'], title=f'Cell Type (resolution={resolution})', save=f"_celltype-{dataset}-{genelist}-{cluster_method}-{resolution}-{embedding_method}.pdf")
                        
                        # 绘制桑基图并保存
                        plot_sankey(adata, cluster_column= cluster_method, celltype_column='celltype', save_path=f"/home/foundation/program/foundation/plot/sankey/sankey-{dataset}-{genelist}-{cluster_method}-{resolution}-{embedding_method}.pdf")
