#! /usr/bin/env python3
# _*_ coding: utf-8 _*_

#@Author:  LiuXing liuxing2@genomics.cn
#@Date:  2023-07-13 14:30:56
#@Last Modified by: LiuXing
#@Last Modified time: 2023-07-13 14:30:56

import os, sys
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import seaborn as sns
import matplotlib.pyplot as plt

def umap(adata: anndata,
         size: int = 10,
         color: str = 'cell2loc_anno',
         legend_loc: str = 'on data',
         legend_fontsize: int = 20,
         figsize: tuple = (6, 6)
         ):
    if 'umap' not in adata.uns['CN']:
        with plt.rc_context({'axes.facecolor':'white','figure.figsize': figsize}):
            sc.pl.umap(adata, 
                       color=[color],
                       size=size, 
                       ncols = 2, 
                       legend_loc=legend_loc, 
                       legend_fontsize=legend_fontsize)
    else:
        if f"{color}_colors" in adata.uns.keys():
            palette = adata.uns[f'{color}_colors']
        else:
            palette = 'tab20'
        umap_result = adata.uns['CN']['umap']
        leiden_cluster = adata.uns['CN']['leiden_cluster']
        fig, ax = plt.subplots(figsize=figsize)
        sns.scatterplot(x = umap_result[:, 0], y = umap_result[:, 1], hue = leiden_cluster, palette = 'tab20', s = 10, ax=ax)
        ax.set_xlabel('UMAP1')
        ax.set_ylabel('UMAP2')
        ax.set_title('Cellular_Neighborhood')
        plt.show()

def spatial(adata: anndata,
            spot_size: int = 100,
            figsize: list = [6, 6]
            ):
    with plt.rc_context({'axes.facecolor':'white','figure.figsize': figsize}):
            sc.pl.spatial(adata, 
                          color=['cell_neighbor'],
                          size=1.3, 
                          img_key='hires', 
                          alpha=1, 
                          spot_size=spot_size)

def heatmap(adata: anndata,
            cmap: str = 'RdYlBu_r',
            figsize: tuple = (5, 3.5),
            row_cluster: bool = False,
            col_cluster: bool = False,
            z_score = None,
            ):
     cn_pct = adata.uns['CN']['cell_composition']
     sns.clustermap(cn_pct,
                    cmap=cmap,
                    figsize = figsize,
                    row_cluster=row_cluster,
                    col_cluster=col_cluster,
                    z_score=z_score)