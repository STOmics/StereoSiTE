#! /usr/bin/env python3
# _*_ coding: utf-8 _*_

#@Author:  LiuXing liuxing2@genomics.cn
#@Date:  2023-07-12 10:34:41
#@Last Modified by: LiuXing
#@Last Modified time: 2023-07-12 10:34:41

import pandas as pd
import os, sys
import anndata

def m2h_homologene(adata: anndata):
    """
    homologous transition gene name from mouse to human

    Parameters
    ----------
    adata
        anndata
    """
    # Biomart tables
    biomart_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "datasets/biomart")
    #h2m_tab = pd.read_csv(os.path.join(biomart_path, "human_to_mouse_biomart_export.csv")).set_index('Gene name')
    m2h_tab = pd.read_csv(os.path.join(biomart_path, "mouse_to_human_biomart_export.csv")).set_index('Gene name')

    hdict = m2h_tab[~m2h_tab['Human gene name'].isnull()]['Human gene name'].to_dict()
    adata.var['original_gene_symbol'] = adata.var_names
    adata.var.index = [hdict[x] if x in hdict.keys() else x.upper() for x in adata.var_names]
    adata.var_names_make_unique()
    return adata

