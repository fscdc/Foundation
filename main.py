import os
os.environ['PYTORCH_CUDA_ALLOC_CONF'] = 'max_split_size_mb:128'

from pathlib import Path
import numpy as np
import scanpy as sc
import warnings
warnings.filterwarnings("ignore", category=Warning)

import scanpy as sc
import scvi
from tools.eval_utils import evaluate_embedding
from tools.tools import load_gene_list, convert_X_to_dense, main_gene_selection, BasicProprecess

from tools.extract_batch import extract_batch
from tools.extract_batch2 import extract_batch2
import argparse
import random

# fix random seed
random.seed(2024)
np.random.seed(2024)

parser = argparse.ArgumentParser(description='MT foundation model')

parser.add_argument('--gene_list', type=str, default='mt genes 1 (1683)', help='which gene list')
parser.add_argument('--embedding_method', type=str, default='PCA', help='which embedding method')
parser.add_argument('--dataset', type=str, default='GSE206785' ,help='which dataset')
parser.add_argument('--batch_list', type=str, default='Patient,Tissue,Platform', help='the obs name which decide an unqiue batch')
parser.add_argument('--cluster_method', type=str, default='leiden', help='which cluster method')
parser.add_argument('--resolution', type=float, default=1.0, help='which resolution')

args = parser.parse_args()

gene_list = []

if args.dataset == 'GSE206785':
    adata_raw = sc.read("/home/foundation/data/foundation/data/data/GSE206785/GSE206785_dense.h5ad")
    adata_raw.obs.rename(columns={"Type": "celltype"}, inplace=True)
    print("After batched:")
    adata_raw = extract_batch(adata_raw, args.batch_list)
    print(adata_raw)
elif args.dataset == 'GSE206785_tumor':
    adata_raw = sc.read("/home/foundation/data/foundation/data/data/GSE206785/GSE206785_dense.h5ad")
    adata_raw.obs.rename(columns={"Type": "celltype"}, inplace=True)
    print("After batched:")
    adata_raw = extract_batch2(adata_raw, args.batch_list)
    print(adata_raw)
elif args.dataset == 'kidney':
    # the h5ad is a batch of kidney data, we did it before (pick one tissue -> convert -> normalize-10000+log -> filter cells)
    adata_raw = sc.read("/home/foundation/data/foundation/data/data/kidney_data_converted2BasicProprecess.h5ad")
    adata_raw.obs.rename(columns={"cell_type": "celltype"}, inplace=True)
    print("Raw:")
    print(adata_raw)
elif args.dataset == 'pancreas':
    # the h5ad is a batch of pancreas data, we did it before (pick one tissue -> convert -> normalize-10000+log -> filter cells -> filter unknown type cells)
    adata_raw = sc.read("/home/foundation/data/foundation/data/data/pancreas_data_converted2BasicProprecess.h5ad")
    adata_raw.obs.rename(columns={"cell_type": "celltype"}, inplace=True)
    print("Raw:")
    print(adata_raw)
else:
    ValueError("Invalid dataset!")
# elif args.dataset == 'GSE261157':
#     adata_raw = sc.read("/home/foundation/data/foundation/data/data/GSE261157/Ctx_adata_dense.h5ad")
#     #adata_raw = sc.read("/home/chenshengquan/data/fengsicheng/Mitochondria/data/GSE261157/Ctx_adata_dense.h5ad")
#     adata_raw.obs.rename(columns={"clusters_annot": "celltype"}, inplace=True)
#     print("Raw:")
#     print(adata_raw)
#     print("After batched:")
#     adata_raw = extract_batch(adata_raw,args.batch_list)
#     print(adata_raw)

adata_raw = BasicProprecess(adata_raw,qc_min_genes=200)
    
import json

if args.gene_list == 'mt_genes':
    with open('/home/foundation/data/foundation/data/data/vocab-1162.json', 'r') as file:
        data_dict = json.load(file)
    gene_list = list(data_dict.keys())
    print("We use 1162 MT genes")
elif args.gene_list == 'all_genes':
    gene_list = adata_raw.var_names
    print("We use all genes")
else:
    ValueError("Invalid gene list!")

adata_new, to_fill_columns = main_gene_selection(adata_raw, gene_list)

print("The number of not found MT genes (Notice: not precise): ",len(to_fill_columns))


if args.gene_list == 'mt_genes':
    pass
else:
    adata_new = adata_raw.copy()


if args.embedding_method == "PCA":
    sc.pp.highly_variable_genes(adata_new, n_top_genes=3000, inplace=True)
    adata_new = adata_new[:, adata_new.var['highly_variable']]
    sc.tl.pca(adata_new, svd_solver="arpack")

    print("PCA: (before evaluation)")
    print(adata_new)

    evaluate_embedding(
        adata_new, 
        embedding_key="X_pca", 
        cluster_method=args.cluster_method, 
        resolution=args.resolution, 
        true_label_key='celltype', 
        args=args,
    )

elif args.embedding_method == "scVI":
    sc.pp.highly_variable_genes(adata_new, n_top_genes=3000, inplace=True)
    adata_new = adata_new[:, adata_new.var['highly_variable']]
    adata_new = adata_new.copy()

    print("scVI: (before evaluation)")
    print(adata_new)

    batch_list = args.batch_list.split(',')
    scvi.model.SCVI.setup_anndata(
        adata_new,
        layer=None,  
        categorical_covariate_keys=batch_list,
    )
    model = scvi.model.SCVI(adata_new, n_layers=2, n_latent=30, gene_likelihood="nb")

    model.train()

    SCVI_LATENT_KEY = "X_scVI"
    adata_new.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()

    evaluate_embedding(
        adata_new, 
        embedding_key="X_scVI", 
        cluster_method=args.cluster_method, 
        resolution=args.resolution, 
        true_label_key='celltype', 
        args=args,
    )
elif args.embedding_method == "scGPT":
    sc.pp.highly_variable_genes(adata_new, n_top_genes=3000, inplace=True)
    adata_new = adata_new[:, adata_new.var['highly_variable']]
    adata_new.var["gene_name"] = adata_new.var_names
    
    model_dir = Path("/home/temporary/data/fengsicheng/scBackdoor/model/scGPT_human")
    import scgpt as scg
    ref_embed_adata = scg.tasks.embed_data(
        adata_new,
        model_dir,
        gene_col="gene_name",
        batch_size=64,
    )
    evaluate_embedding(
        ref_embed_adata, 
        embedding_key="X_scGPT", 
        cluster_method=args.cluster_method, 
        resolution=args.resolution, 
        true_label_key='celltype', 
        args=args,
    )

elif args.embedding_method == "scGPT_allgene":
    adata_new.var["gene_name"] = adata_new.var_names
    
    model_dir = Path("/home/temporary/data/fengsicheng/scBackdoor/model/scGPT_human")
    import scgpt as scg
    ref_embed_adata = scg.tasks.embed_data(
        adata_new,
        model_dir,
        gene_col="gene_name",
        batch_size=64,
    )
    evaluate_embedding(
        ref_embed_adata, 
        embedding_key="X_scGPT", 
        cluster_method=args.cluster_method, 
        resolution=args.resolution, 
        true_label_key='celltype', 
        args=args,
    )

print("Done for calculating!")
print("Start to draw plot!")