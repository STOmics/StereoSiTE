#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

"""
@File        :   tf_infer.py
@Description :   tf_infer.py
@Author      :   Liuchuandong liuchuandong@genomics.cn
@Date        :   2023/08/07 10:26:24
"""
import decoupler as dc
import pandas as pd

def ulm(de_result,grn,weight,save):
    '''
    Inferring transcription factor (TF) activities from our transcriptomics data by using decoupler's Univariate Linear Model(ulm) method.
    
    Parameters
    ----------
    de_result
        dataframe of deseq2 result.
    grn
        Gene regulatory network from CollecTRI.
    weight
        Column name in GRN net with weights.
    save
        Path to where to save the plot. Infer the filetype if ending on {`.pdf`, `.png`, `.svg`}.
    ----------

    Returns
        dataframe contains TF enrichment scores and pvalues
    '''
    
    mat = de_result[['stat']].T
    tf_acts, tf_pvals = dc.run_ulm(mat=mat, net=grn, weight=weight,verbose=True)
    tf_df = pd.DataFrame({'tf_acts':tf_acts.T['stat'],
                          'tf_pvals':tf_pvals.T['stat']},index=tf_acts.columns)
    dc.plot_barplot(tf_acts, 'stat', top=25, vertical=True,save=save)

    return(tf_df)