#! /usr/bin/env python3
# _*_ coding: utf-8 _*_

#@Author:  LiuXing liuxing2@genomics.cn
#@Date:  2023-07-12 15:29:29
#@Last Modified by: LiuXing
#@Last Modified time: 2023-07-12 15:29:29

import os, sys
import anndata
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
from scipy import sparse
import numpy as np
import pandas as pd
from tqdm.notebook import tqdm

def _intensity_show(LRadata: anndata,
                    cells: list,
                    genes: list,
                    l: np.ndarray,
                    r: np.ndarray,
                    key: str = 'spatial_exp'):
    """
    draw the spatial cell interaction intensity in space.
    """
    df = pd.DataFrame()
    df['x'] = LRadata.obsm["spatial"][:,0].astype(int)
    df['y'] = LRadata.obsm["spatial"][:,1].astype(int)
    df['intensity'] = LRadata.obsp[key].sum(axis=1).A1*(LRadata.obs['cell2loc_anno'].isin(cells).any())
    df[genes[0]] = l
    df[genes[1]] = r
    colors = ['tab:blue', 'tab:orange', 'tab:red']
    fig, ax = plt.subplots(figsize=(4,5), dpi=200)
    i=0
    for gene in genes:
        df1 = df[df[gene]>0]
        scatter = ax.scatter(x =  df1['x'], y = df1['y'], c = colors[i], label = gene, edgecolors='none', alpha = 0.5, s = 1)
        i+=1
    if (df.shape[0] > 0):
        scatter = ax.scatter(x = df['x'], y = df['y'], c = colors[2], s = df['intensity'], edgecolors='none', label='intensity', alpha = 0.5)
        handles, labels = scatter.legend_elements(prop="sizes", alpha=0.6)
        legend1 = ax.legend(handles, labels, loc = 'upper left', prop = {'size': 3}, title = 'intensity')
        ax.add_artist(legend1)
    else:
        return 0
    ax.set_ylim(bottom = df['y'].max(), top = 0)
    ax.legend(prop = {'size': 6})
    plt.gca().set_aspect(1)
    plt.axis('off')
    plt.title(f"{cells[0]}_{cells[1]} ({genes[0]}_{genes[1]})")
    
def intensity_insitu(adata: anndata,
                  cells: list,
                  genes: list,
                  anno: str= 'cell2loc_anno',
                  radius: int = None,
                  connectivities_key: str = "spatial_connectivities",
                  complex_process_model = 'mean') -> int:
    """
    Calculate the spatial cell interaction intensity between specified cells and ligand receptor genes.

    Parameters
    ----------
    adata
        anndata
    cells
        list contains sender and receiver cell type. [sender_cell, receiver_cell]
    genes
        list contains ligand and receptor genes. [ligand, receptor]
    anno
        annotation key, default=cell2loc_anno.
    radius
        radius threshold when constructing nearest neighbor graph. Only be used when the neighbor graph doesn't exist
    connectivities_key
        obtain the constructed nearest neighbor graph in adata.obsp[connectivities_key]. If this doesn't exist, construct
        the neighbor graph with specified radius
    complex_process_model
        determine how to deal with the complexed ligand and receptor which contain multiple subunits. There are two options: mean, min.
        mean: calculate the mean expression of all subunits to represent the complex
        min: pick the minimal expression of all subunits to represent the complex
    
    Returns
    ----------
        intensity value of the entire slide
    """ 
    exp_key = "spatial_exp"
    if "_" in genes[0]:
        if complex_process_model == 'min':
            l = (adata.obs[anno]==cells[0])*(adata[:,genes[0].split("_")].X.min(axis=1).toarray()[:,0])
        elif complex_process_model == 'mean':
            l = (adata.obs[anno]==cells[0])*(adata[:,genes[0].split("_")].X.mean(axis=1).A1)
        else:
            raise Exception(f"complex_process_model should be mean or min, but got {complex_process_model}.")
    else:
        l = (adata.obs[anno]==cells[0])*(adata[:,genes[0]].X.sum(axis=1).A1)
    if "_" in genes[1]:
        if complex_process_model == 'min':
            r = (adata.obs[anno]==cells[1])*(adata[:,genes[1].split("_")].X.min(axis=1).toarray()[:,0])
        elif complex_process_model == 'mean':
            r = (adata.obs['cell2loc_anno']==cells[1])*(adata[:,genes[1].split("_")].X.mean(axis=1).A1)
        else:
            raise Exception(f"complex_process_model should be mean or min, but got {complex_process_model}.")
    else:
        r = (adata.obs[anno]==cells[1])*(adata[:,genes[1]].X.sum(axis=1).A1)
    
    if connectivities_key in adata.obsp.keys() and radius == None:
        connect_matrix = adata.obsp[connectivities_key]
    elif connectivities_key not in adata.obsp.keys() and radius > 10 :
        key_added = connectivities_key.replace("_connectivities", "")
        sq.gr.spatial_neighbors(adata, radius=radius*2, coord_type="generic", key_added=key_added)
        connect_matrix = adata.obsp[connectivities_key]
    else:
        raise Exception(f"The connectivities_key ({connectivities_key}) dosn't exist in adata.obsp, and radius has not be specified with a value >= 10")

    l = l.values
    r = r.values
    l_rows = np.where(l > 0)[0]
    r_cols = np.where(r > 0)[0]
    dst = np.where(connect_matrix[l_rows,:][:,r_cols].todense()>0)
    exps = l[l_rows[dst[0]]] + r[r_cols[dst[1]]]
    spatial_exp = sparse.csr_matrix((exps, (l_rows[dst[0]], r_cols[dst[1]])), shape=connect_matrix.shape, dtype=int)
    adata.obsp[exp_key] = spatial_exp
    _intensity_show(adata, cells, genes, l, r, key=exp_key)
    return spatial_exp.sum()

def intensities_with_radius(adata, pairs = None):
    """
    draw the line plot with radius as x and intensity per area as y

    Parameters
    ---------------
    adata
        anndata
    pairs
        pair list that will be drawed in the plot
        for example Teff|Microphage(APP|CD74)
    """
    try:
        plot_df = adata.uns['intensities_with_radius']
    except Exception:
        print("there must be intensities_with_radius in adata.uns, please run intensity.intensities_with_radius before draw this plot.")
    if pairs == None:
        columns = plot_df.columns.values.tolist()
        columns.remove('radius')
    else:
        columns = pairs

    import matplotlib.colors as mcolors
    colors = list(mcolors.TABLEAU_COLORS.keys())
    i = 0
    for column in columns:
        x = plot_df['radius']
        y = plot_df[column]/(3.14*(plot_df['radius']**2))
        plt.plot(x, y,
                 linestyle="-", linewidth=2, color=colors[i],
                 marker='o', markersize=4, markeredgecolor='black', markerfacecolor=colors[i],
                 label=column,
                 )
        i+=1
    plt.xlabel('radius (µm)')
    plt.ylabel('intensity / area')
    plt.legend()

