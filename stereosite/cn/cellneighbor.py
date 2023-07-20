#! /usr/bin/env python3
# _*_ coding: utf-8 _*_

#@Author:  LiuXing liuxing2@genomics.cn
#@Date:  2023-07-12 15:24:52
#@Last Modified by: LiuXing
#@Last Modified time: 2023-07-12 15:24:52



import os, sys
from optparse import OptionParser
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import sparse as sp
import logging as logg

def cn_deconvolve(adata: anndata,
                     use_rep: str = 'q05_cell_abundance_w_sf',
                     n_neighbors: int = 20,
                     resolution: float = 0.4,
                     min_dist: float = 0.2,
                     random_stat: int = 100,
                     key_added: str = None):
    """
    Cluster bins with similar cellular composition, which was calculated from cell type deconvolution result.

    Parameters
    ----------
    adata
        anndata
    use_rep
        Key value of anndata.obsm, can be used to obtain deconvolution matrix of each bin
    n_neighbors
        Number of neighbors that will be connected with the center point when constructing k-nearest neighbor graph
    resolution
        A parameter value controlling the coarseness of the clustering. Higher values lead to more clusters. Used for leiden function.
    min_dist
        The effective minimum distance between embedded points. Smaller values will result in a more clustered/clumped embedding where 
        nearby points on the manifold are drawn closer together, while larger values will result on a more even dispersal of points. 
        The value should be set relative to the spread value, which determines the scale at which embedded points will be spread out. 
        The default of in the umap-learn package is 0.1.
    key_added
        Column name of adata.obs, which store the result of leiden clustering. 
    """
    if key_added == None:
        key_added = 'cell_neighbor'
    sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=n_neighbors, key_added='CN')
    sc.tl.leiden(adata, resolution=resolution, key_added=key_added, neighbors_key='CN')
    sc.tl.umap(adata, min_dist=min_dist, spread=1, random_state=random_stat, neighbors_key='CN')

    # analysis each Cellular Neighborhood's celltype composition according to bin level cell2location annotation matrix
    celltypes= [str[23:] for str in list(adata.obsm[use_rep])]
    # sum celltypes in each Cellular Neighborhood
    cn_pct = pd.pivot_table(adata.obs,columns=key_added,values=celltypes,aggfunc=np.sum).T
    # calculate percentage
    cn_pct = cn_pct.apply(lambda x: x/x.sum(), axis=1)
    adata.uns['CN']['cell_composition'] = cn_pct
    
def cn_cellbin(adata: anndata,
               bin_size: int = 200,
               anno: str = 'cell2loc_anno',
               n_neighbors: int = 20,
               resolution: float = 0.4,
               min_dist: float = 0.1,
               key_added: str = None,
               random_stat:int = 100,
               ):
    """
    Bin the annotated ST data in single-resolution, and cluster bins with similar cellular composition.

    Parameters
    ----------
    adata
        anndata
    bin_size
        Determines the size of bin. Take 200 as an example, calculate the cellular composition in every sparately arranged 200x200 square bin.
    anno
        Column name of anndata.obs, which contains the cell type information.
    n_neighbors
        Number of neighbors that will be connected with the center point when constructing k-nearest neighbor graph
    resolution
        A parameter value controlling the coarseness of the clustering. Higher values lead to more clusters. Used for leiden function.
    min_dist
        The effective minimum distance between embedded points. Smaller values will result in a more clustered/clumped embedding where 
        nearby points on the manifold are drawn closer together, while larger values will result on a more even dispersal of points. 
        The value should be set relative to the spread value, which determines the scale at which embedded points will be spread out. 
        The default of in the umap-learn package is 0.1.
    key_added
        Column name of adata.obs, which store the result of leiden clustering. 
    """
    if key_added == None:
        key_added = 'cell_neighbor'
    bin_cor =[str(x[0]) + "-" + str(x[1]) for x in ((adata.obsm['spatial']//bin_size)*bin_size+(bin_size/2)).astype(int)]
    adata.obs['bin_cor'] = bin_cor
    groups = adata.obs[anno].groupby(adata.obs['bin_cor'])
    cellbin_count = pd.DataFrame(columns = adata.obs['cell2loc_anno'].unique())
    for group in groups:
        pct = group[1].value_counts()
        cellbin_count.loc[group[0]] = pct.values
    
    from sklearn.neighbors import NearestNeighbors
    #construct nearest neighbor graph
    nbrs = NearestNeighbors(n_neighbors=n_neighbors, algorithm='ball_tree').fit(cellbin_count)
    distances, indices = nbrs.kneighbors(cellbin_count)
    X = sp.coo_matrix(([], ([], [])), shape=(indices.shape[0], 1))
    adjacency = _get_sparse_matrix_from_indices_distances_umap(indices, distances, distances.shape[0], n_neighbors)
    g = _get_igraph_from_adjacency(adjacency, directed=False)

    #cluster bins by leiden
    import leidenalg
    partition_type = leidenalg.RBConfigurationVertexPartition
    partition_kwargs = dict()
    partition_kwargs['weights'] = np.array(g.es['weight']).astype(np.float64)
    partition_kwargs['n_iterations'] = -1
    partition_kwargs['seed'] = 1
    partition_kwargs['resolution_parameter'] = resolution
    part = leidenalg.find_partition(g, partition_type, **partition_kwargs)
    groups = np.array(part.membership)
    cn_cluster = dict(zip(cellbin_count.index.values, groups))
    adata.obs[key_added] = adata.obs['bin_cor'].map(cn_cluster).astype('category')
    adata.uns['CN'] = dict()
    adata.uns['CN']['params'] = {'n_neighbors': n_neighbors,
                                 'random_state': random_stat,
                                 }
    adata.uns['CN']['leiden_cluster'] = groups
    values = cellbin_count.columns
    cellbin_count[key_added] = groups
    cn_pct = pd.pivot_table(cellbin_count, columns = key_added, values = values, aggfunc=np.sum).T
    cn_pct = cn_pct.apply(lambda x: x/x.sum()*100, axis=1)
    adata.uns['CN']['cell_composition'] = cn_pct

    #umap
    import umap
    umap_model = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, n_components=2, random_state=random_stat)
    umap_result = umap_model.fit_transform(cellbin_count.values)
    
    adata.uns['CN']['umap'] = umap_result


def _get_sparse_matrix_from_indices_distances_umap(
    knn_indices, knn_dists, n_obs, n_neighbors
):
    rows = np.zeros((n_obs * n_neighbors), dtype=np.int64)
    cols = np.zeros((n_obs * n_neighbors), dtype=np.int64)
    vals = np.zeros((n_obs * n_neighbors), dtype=np.float64)

    for i in range(knn_indices.shape[0]):
        for j in range(n_neighbors):
            if knn_indices[i, j] == -1:
                continue  # We didn't get the full knn for i
            if knn_indices[i, j] == i:
                val = 0.0
            else:
                val = knn_dists[i, j]

            rows[i * n_neighbors + j] = i
            cols[i * n_neighbors + j] = knn_indices[i, j]
            vals[i * n_neighbors + j] = val

    result = sp.coo_matrix((vals, (rows, cols)), shape=(n_obs, n_obs))
    result.eliminate_zeros()
    return result.tocsr()

def _get_igraph_from_adjacency(adjacency, directed=None):
    """Get igraph graph from adjacency matrix."""
    import igraph as ig

    sources, targets = adjacency.nonzero()
    weights = adjacency[sources, targets]
    if isinstance(weights, np.matrix):
        weights = weights.A1
    g = ig.Graph(directed=directed)
    g.add_vertices(adjacency.shape[0])  # this adds adjacency.shape[0] vertices
    g.add_edges(list(zip(sources, targets)))
    try:
        g.es['weight'] = weights
    except KeyError:
        pass
    if g.vcount() != adjacency.shape[0]:
        logg.warning(
            f'The constructed graph has only {g.vcount()} nodes. '
            'Your adjacency matrix contained redundant nodes.'
        )
    return g
          
