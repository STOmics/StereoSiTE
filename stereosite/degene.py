#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

"""
@File        :   de_gene.py
@Description :   de_gene.py
@Author      :   Liuchuandong liuchuandong@genomics.cn
@Date        :   2023/08/04 17:08:24
"""
import decoupler as dc
import scanpy as sc
import anndata
import sys, os
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# Only needed for visualization:
import matplotlib.pyplot as plt
import seaborn as sns

def deseq2(adata,sample_col='sample',groups_col='cell_neighbor',
           contrast=['5','Others'],batch=None,save=None):
    '''
    Perform differential gene expression analysis in pseudo-bulk manner by using 'decoupler' and 'pydeseq2' modules.
    Here just give a one line method in StereoSite workflow.

    Parameters
    ----------
    adata
        Anndata objects, make sure .X is raw counts.
    sample_col
        Column of obs contains samples names.
    groups_col
        Column of obs contains groups names.
    contrast
        list contains group names you want to compare,['5','Others'] means CN5 vs Others_CN, while ['5','1'] means CN5 vs CN1
    batch
        Column of obs contains batch names 
    save
        Path to where to save the volcano plot. Infer the filetype if ending on {`.pdf`, `.png`, `.svg`}.
    ----------
    
    Returns
        dataframe of deseq2 result
    
    '''
    # get pseudobulk
    pdata = dc.get_pseudobulk(
        adata,
        sample_col=sample_col,
        groups_col=groups_col,
        mode='sum',
        min_cells=0,
        min_counts=0)
    # Filter genes
    genes = dc.filter_by_expr(pdata, group=groups_col, min_count=0, min_total_count=10)
    pdata = pdata[:, genes].copy()
    if 'Others' in contrast:
        pdata.obs['deGroup'] = [contrast[0] if i==contrast[0] else 'Others' for i in pdata.obs[groups_col]]
    else:
        pdata.obs['deGroup'] = pdata.obs[groups_col]
    # Build DESeq2 object
    if batch is None:
        design_factors = ['deGroup']
    else:
        design_factors=['deGroup',batch]
    dds = DeseqDataSet(
        adata=pdata,
        design_factors=design_factors,
        ref_level=['deGroup', contrast[1]],
        refit_cooks=True,
        n_cpus=8,
    )
    # Compute LFCs
    dds.deseq2()
    # Extract contrast between CN5 vs Others
    stat_res = DeseqStats(dds, contrast=["deGroup", contrast[0], contrast[1]], n_cpus=8)
    # Compute Wald test
    stat_res.summary()
    # Shrink LFCs
    stat_res.lfc_shrink(coeff='deGroup_'+contrast[0]+'_vs_'+contrast[1])
    # Extract results
    results_df = stat_res.results_df
    dc.plot_volcano_df(results_df, x='log2FoldChange', y='padj', top=20,save=save)
    return(results_df)
