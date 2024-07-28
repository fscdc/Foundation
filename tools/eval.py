import numpy as np
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.metrics.cluster import fowlkes_mallows_score
from sklearn.metrics.cluster import homogeneity_score
from sklearn.metrics.cluster import adjusted_mutual_info_score
from sklearn.metrics.cluster import completeness_score
import sklearn
import sklearn.neighbors
import scanpy as sc
from sklearn.model_selection import cross_validate
from sklearn.svm import SVC
from sklearn.metrics import cohen_kappa_score, make_scorer

def svm_cross_validation(mtx, target, Kfold=5):
    """
    Perform SVM cross validation for the given matrix and target.

    Parameters
    ----------
    mtx
        Matrix of shape (n_samples, n_features)
    target
        Target labels of shape (n_samples,)
    Kfold
        Number of folds for cross validation
    """
    target = np.array(target)
    target_unique = np.unique(target).reshape(1, -1)
    target_onehot = (target.reshape(-1, 1)==target_unique).astype(int)
    target = target_onehot.argmax(-1)
    svc = SVC()
    cv_results = cross_validate(svc, mtx, target,
                                scoring=("accuracy", "f1_macro", "f1_weighted"),
                                cv=Kfold, n_jobs=Kfold)
    svc = SVC()
    kappa_score = make_scorer(cohen_kappa_score)
    kappa = cross_validate(svc, mtx, target,
                                scoring=kappa_score,
                                cv=Kfold, n_jobs=Kfold)["test_score"]
    
    acc, kappa, mf1, wf1 = cv_results["test_accuracy"].mean(), kappa.mean(), cv_results["test_f1_macro"].mean(), cv_results["test_f1_weighted"].mean()
    
    print('Accuracy: %.3f, Kappa: %.3f, mF1: %.3f, wF1: %.3f' % (acc, kappa, mf1, wf1))
    
    return acc, kappa, mf1, wf1

def cluster_metrics(target, pred):
    """
    calculate ARI, AMI, NMI, FMI, Completeness, Homogeneity

    Parameters
    ----------
    target
        True labels
    pred
        Predicted labels
    """
    target = np.array(target)
    pred = np.array(pred)
    
    ari = adjusted_rand_score(target, pred)
    ami = adjusted_mutual_info_score(target, pred)
    nmi = normalized_mutual_info_score(target, pred)
    fmi = fowlkes_mallows_score(target, pred)
    comp = completeness_score(target, pred)
    homo = homogeneity_score(target, pred)
    print('ARI: %.3f, AMI: %.3f, NMI: %.3f, FMI: %.3f, Comp: %.3f, Homo: %.3f' % (ari, ami, nmi, fmi, comp, homo))
    
    return ari, ami, nmi, fmi, comp, homo

def mean_average_precision(x: np.ndarray, y: np.ndarray, k: int=30, **kwargs) -> float:
    r"""
    Mean average precision
    Parameters
    ----------
    x
        Coordinates
    y
        Cell_type/Layer labels
    k
        k neighbors
    **kwargs
        Additional keyword arguments are passed to
        :class:`sklearn.neighbors.NearestNeighbors`
    Returns
    -------
    map
        Mean average precision
    """
    
    def _average_precision(match: np.ndarray) -> float:
        if np.any(match):
            cummean = np.cumsum(match) / (np.arange(match.size) + 1)
            return cummean[match].mean().item()
        return 0.0
    
    y = np.array(y)
    knn = sklearn.neighbors.NearestNeighbors(n_neighbors=min(y.shape[0], k + 1), **kwargs).fit(x)
    nni = knn.kneighbors(x, return_distance=False)
    match = np.equal(y[nni[:, 1:]], np.expand_dims(y, 1))
    
    return np.apply_along_axis(_average_precision, 1, match).mean().item()

def rep_metrics(adata, use_rep, key, k_map=30):
    """
    Calculate MAP, ASW, and cLISI for a given embedding.

    Parameters
    ----------
    adata
        AnnData object
    use_rep
        Key to access the embedding in adata.obsm
    key
        Key to access the true labels in adata.obs
    k_map
        Number of neighbors for MAP calculation
    """
    import scib
    if key not in adata.obs or use_rep not in adata.obsm:
        print("KeyError")
        return None
    
    adata.obs[key] = adata.obs[key].astype("category")
    MAP = mean_average_precision(adata.obsm[use_rep].copy(), adata.obs[key], k=k_map)
    ASW = scib.me.silhouette(adata, label_key=key, embed=use_rep)
    cLISI = scib.me.clisi_graph(adata, label_key=key, type_="embed", use_rep=use_rep)
    print('MAP: %.3f, ASW: %.3f, cLISI: %.3f' % (MAP, ASW, cLISI))
    
    return MAP, ASW, cLISI

def get_N_clusters(adata, n_cluster, cluster_method='louvain', range_min=0, range_max=3, max_steps=30, tolerance=0,usedefaultres=None,resolution=1.0):
    """
    Tune the resolution parameter in clustering to make the number of clusters and the specified number as close as possible.
    
    Parameters
    ----------
    adata
        AnnData object of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.
    n_cluster
        Specified number of clusters.
    cluster_method
        Method (`louvain` or `leiden`) used for clustering. By default, cluster_method='louvain'.
    range_min
        Minimum clustering resolution for the binary search.
    range_max
        Maximum clustering resolution for the binary search.
    max_steps
        Maximum number of steps for the binary search.
    tolerance
        Tolerance of the difference between the number of clusters and the specified number.

    Returns
    -------
    adata
        AnnData object with clustering assignments in `adata.obs`:

        - `adata.obs['louvain']` - Louvain clustering assignments if `cluster_method='louvain'`.
        - `adata.obs['leiden']` - Leiden clustering assignments if `cluster_method='leiden'`.

    """ 
    if usedefaultres==False:
        this_step = 0
        this_min = float(range_min)
        this_max = float(range_max)
        while this_step < max_steps:
            this_resolution = this_min + ((this_max-this_min)/2)
            if cluster_method=='leiden':
                sc.tl.leiden(adata, resolution=this_resolution)
            if cluster_method=='louvain':
                sc.tl.louvain(adata, resolution=this_resolution)
            this_clusters = adata.obs[cluster_method].nunique()

            if this_clusters > n_cluster+tolerance:
                this_max = this_resolution
            elif this_clusters < n_cluster-tolerance:
                this_min = this_resolution
            else:
                print("Succeed to find %d clusters at resolution %.3f."%(n_cluster, this_resolution))
                return adata
            this_step += 1

        print('Cannot find the number of clusters.')
        return adata
    else:
        this_resolution=resolution
        if cluster_method=='leiden':
            sc.tl.leiden(adata, resolution=this_resolution)
        if cluster_method=='louvain':
            sc.tl.louvain(adata, resolution=this_resolution)
        return adata

def save_adata(adata, save_path, file_name):
    """
    Save an AnnData object to a file.

    Parameters:
    - adata: AnnData object
    - save_path: path to save the file
    - file_name: name of the file
    """
    sc.write(save_path + file_name, adata)
    print(f"Saved adata to {save_path + file_name}")

def evaluate_embedding(adata, embedding_key=None, cluster_method='leiden', resolution=1.0, true_label_key='celltype', args=None):
    """
    Evaluate embedding method by clustering and calculating ARI and AMI.

    Parameters:
    - adata: AnnData object
    - embedding_key: key to access the embedding in adata.obsm
    - cluster_method: clustering method ('leiden' or 'louvain')
    - resolution: resolution parameter for clustering
    - true_label_key: key of true labels in adata.obs

    Returns:
    - result_dict: dictionary with ARI and AMI scores
    """
    
    # if resolution is 888.8, use binary search to find the best resolution
    # else use the given resolution
    usedefaultres=True
    if resolution==888.8:
        usedefaultres=False

    # Perform clustering
    sc.pp.neighbors(adata, use_rep=embedding_key)
    real_cluster = len(adata.obs['celltype'].values.unique())
    
    if cluster_method == 'leiden':
        adata=get_N_clusters(adata, n_cluster=real_cluster, cluster_method='leiden', range_min=0, range_max=3, max_steps=30, tolerance=0, usedefaultres=usedefaultres,resolution=resolution)
        cluster_labels = adata.obs['leiden']
    elif cluster_method == 'louvain':
        adata=get_N_clusters(adata, n_cluster=real_cluster, cluster_method='louvain', range_min=0, range_max=3, max_steps=30, tolerance=0, usedefaultres=usedefaultres,resolution=resolution)
        cluster_labels = adata.obs['louvain']
    else:
        raise ValueError("Invalid clustering method. Choose 'leiden' or 'louvain'.")
    
    # True labels
    true_labels = adata.obs[true_label_key]
    
    # Calculate Metrics for cluster
    cluster_metrics(target=true_labels, pred=cluster_labels)

    # save adata
    save_path2 = "/home/foundation/program/foundation-new/record/temp-h5ad/"
    file_name = f"result-{args.dataset}-{args.gene_list}-{args.cluster_method}-{args.resolution}-{args.embedding_method}.h5ad"
    save_adata(adata, save_path2, file_name)

    
