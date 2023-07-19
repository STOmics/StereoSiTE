#! /usr/bin/env python3
# _*_ coding: utf-8 _*_

#@Author:  LiuXing liuxing2@genomics.cn
#@Date:  2023-07-12 11:01:58
#@Last Modified by: LiuXing
#@Last Modified time: 2023-07-12 11:01:58

 

import sys, os

import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import csv
from scipy import sparse, stats
import gzip
import tifffile
import cv2
import tifffile

class Gem_Reader():
    """
    class for transfer gene expression matrix to anndata

    Parameters
    ----------
    gem_file
        gene expression matrix (gem) file path. 
        Gem file should include at least 4 columns: geneID, x, y, MIDCount. MIDCount indicates gene expression count, and can be specified by parameter count_key
        If cell segmentation has been done, gen file is consist of 5 columns: geneID, x, y, MIDCount, label. transcripts with same label belong to the same single cell.
    tissue_mask
        tissue mask file path, should be tiff format, used to extract gene expression data under the tissue covered region
    count_key
        column name of gene expression count in gem file
    cell_label_key
        column name of cell labels in gem file with single-cell segmentation.
    """
    def __init__(self,
                 gem_file: str,
                 tissue_mask: str = None,
                 count_key: str = 'MIDCount',
                 cell_label_key: str = 'label'):
        self.gem_file = gem_file
        self.tissue_mask = tissue_mask
        self.count_key = count_key
        self.cell_label_key = cell_label_key
        self._read_gem()

    def _read_gem(self) -> pd.DataFrame:
        """
        Read gene expression matrix from the given file. And extract data under tissue covered region if tissue mask file was given
        """

        columntypes = {"geneID": 'category', 
                       "x": int, 
                       'y': int, 
                       self.count_key: int,
                       }
    
        gem = pd.read_csv(self.gem_file, sep="\t", quoting=csv.QUOTE_NONE, comment="#", dtype=columntypes)

        gem.rename(columns={self.count_key : 'counts',
                            self.cell_label_key: 'label'})
        
        #extract gene expression under tissue covered region
        if self.tissue_mask != None:
            tissue_mask = tifffile.imread(self.tissue_mask)
            maxY, maxX = tissue_mask.shape
            if maxX > gem['x'].max() or maxY > gem['y'].max():
                print ("WARMING: mask is out of bounds")
                gem = gem.loc[(gem['x'] < maxX)&(gem['y'] < maxY)]
            gem = gem.loc[tissue_mask[gem['y'], gem['x']] > 0]
        self.gem = gem

    def gem2anndata(self, bin_size=50) -> anndata:
        """
        transfer gem to anndata

        Parameters
        ----------
        bin_size
            Specify bin size with this parameter. For example, bin_sie = 50 means bining spots in the same 50X50 square.
        
        Returns
        ----------
        anndata
            contains gene expression vector and spatial coordinate of each bin.
        """
    
        half_bin_size = int(bin_size/2)
    
        self.gem['x'] = (self.gem['x']//bin_size)*bin_size + half_bin_size
        self.gem['y'] = (self.gem['y']//bin_size)*bin_size + half_bin_size 
        self.gem['cell'] = self.gem['x'].astype(str) + "-" + self.gem['y'].astype(str)

        cells = self.gem['cell'].unique()
        genes = self.gem['geneID'].unique()
    
        cells_dict = dict(zip(cells, range(0, len(cells))))
        genes_dict = dict(zip(genes, range(0, len(genes))))
        rows = self.gem['cell'].map(cells_dict)
        cols = self.gem['geneID'].map(genes_dict)
        expMtx = sparse.csr_matrix((self.gem[self.count_key].values, (rows, cols)), shape=(cells.shape[0], genes.shape[0]), dtype=np.int32)

        obs = pd.DataFrame(index = cells)
        var = pd.DataFrame(index = genes)
        adata = anndata.AnnData(X = expMtx, obs = obs, var = var)
        positions = np.array(list(map(lambda x: [int(v) for v in x.strip().split("-")], adata.obs.index)))
        adata.obs['x'] = positions[:,0]
        adata.obs['y'] = positions[:,1]
        adata.obsm['spatial'] = adata.obs[['x', 'y']].values
        return adata

    def cellbin2anndata(self) -> anndata:
        """
        transfer gem with single-cell segmentation to anndata

        Returns
        -----------
        anndata
            contains gene expression vector and sptial coordinate of each cell 
        """

        gem = self.gem[self.gem[self.cell_label_key]!=0]
    
        cells = gem[self.cell_label_key].unique()
        genes = gem['geneID'].unique()
        cells_dict = dict(zip(cells, range(0, len(cells))))
        genes_dict = dict(zip(genes, range(0, len(genes))))
        rows = gem[self.cell_label_key].map(cells_dict)
        cols = gem['geneID'].map(genes_dict)

        expMtx = sparse.csr_matrix((gem[self.count_key].values, (rows, cols)), shape=(cells.shape[0], genes.shape[0]), dtype=np.int32)

        obs = pd.DataFrame(index = cells)
        var = pd.DataFrame(index = genes)
        adata = anndata.AnnData(X = expMtx, obs = obs, var = var)
        spatialgroup = gem[['x', 'y']].groupby(gem[self.cell_label_key])
        spatialdf = spatialgroup.agg(lambda x: (x.max()+x.min())/2)
        spatialdf = spatialdf.reset_index()
        spatialdf[self.cell_label_key] = spatialdf[self.cell_label_key].astype('category')
        spatialdf[self.cell_label_key].cat.reorder_categories(cells, inplace=True)
        spatialdf.sort_values(self.cell_label_key, inplace=True)
        adata.obsm['spatial'] = spatialdf[['x', 'y']].values
        return adata

    def gem_with_cellmask_2anndata(self, cell_mask: str) -> anndata:
        """
        extract single-cell gem based on the cell mask and transfer it to anndata

        Parameters
        -----------
        cell_mask
            cell mask file path, should be tiff format, used to extract gene expression of each single-cell
        
        Returns
        -----------
        anndata
            contains gene expression vector and sptial coordinate of each cell
        """

        mask = cv2.imread(cell_mask, -1)
        _, labels = cv2.connectedComponents(mask)
        tissuedf = pd.DataFrame()
        dst = np.nonzero(labels)

        tissuedf['x'] = dst[1]
        tissuedf['y'] = dst[0]
        tissuedf[self.cell_label_key] = labels[dst]

        res = pd.merge(self.gem, tissuedf, on=['x', 'y'], how='inner')
        self.gem = res
        adata = self.cellbin2anndata()
        return adata
