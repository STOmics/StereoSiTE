#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

"""
@File        :   get_hub_net.py
@Description :   get_hub_net.py
@Author      :   Liuchuandong liuchuandong@genomics.cn
@Date        :   2023/08/01 11:42:43
"""



import os, sys
from optparse import OptionParser
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

POS_SEED=99
HUB_NODE_SIZE = 366
BASIC_NODE_SIZE = 100
PPI_outdir = 'PPI'
os.makedirs(PPI_outdir, exist_ok=True)
TF_outdir = 'TF'
os.makedirs(TF_outdir, exist_ok=True)

def _edge_colors(graph,cmap='binary',min_weight=0.4,max_weight=1,unit=0.001):
    '''
    set edges color by STRING combined_score
    '''
    norm=matplotlib.colors.Normalize(vmin=min_weight, vmax=max_weight)
    scale_map = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
    colors = []
    for i in graph.edges():
        weight = graph.edges[i]['combined_score']
        colors.append( matplotlib.colors.rgb2hex(scale_map.to_rgba(weight*unit)) )
    return({'color_list':colors,'scale_colormap':scale_map})

def _node_colors(graph,logFC_df,name,cmap='Reds',min_weight=None,max_weight=None,unit=0.1):
    '''
    set node colors by logFC
    logFC_df
        Dataframe contains gene's logFoldChange column,and node names as row index
    name
        Column name in logFC_df with gene's logFoldChange
    '''
    logFC_df = logFC_df.loc[list(graph.nodes()),:]
    for i in graph.nodes():
        graph.nodes[i]['logFC']=logFC_df.loc[i,name]
    if min_weight==None:
        min_weight= round(min(logFC_df[name]),1)
    if max_weight==None:
        max_weight= round(max(logFC_df[name]),1)
        
    norm=matplotlib.colors.Normalize(vmin=min_weight, vmax=max_weight)
    scale_map = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
    colors = []
    for i in graph.nodes():
        fc = round(graph.nodes[i]['logFC'],1)
        colors.append( matplotlib.colors.rgb2hex(scale_map.to_rgba(fc)) )
    return({'color_list':colors,'scale_colormap':scale_map})


def ppi_hub_net(hub_net,hub_gene,logFC_df,name,save=PPI_outdir+'hub_net_cluster1.pdf',figsize=(6,6.8),hspace = 0.4,
                 edge_min_weight=0.8,edge_max_weight=1.0,edge_unit=0.001,edge_cmap='binary',
                 node_min_weight=0,node_max_weight=None,node_unit=0.1,node_cmap='Reds'):
    '''
    plot ppi network

    Parameters
    ----------
    hub_net
        Networkx object, the net need to plot.
    hub_gene
        The hub genes of ppi net.
    logFC_df
        Differential analysis result contains gene's logFoldChange column.
    name
        Column name in logFC_df with gene's logFoldChange
    ----------
    '''
    fig, ax=plt.subplots(3,1,figsize=figsize,height_ratios=[28,1,1])
    edge_colors = _edge_colors(hub_net,cmap=edge_cmap,min_weight=edge_min_weight,max_weight=edge_max_weight,unit=edge_unit)
    node_colors = _node_colors(hub_net,logFC_df=logFC_df,name=name,cmap=node_cmap,min_weight=node_min_weight,max_weight=node_max_weight,unit=node_unit)
    pos = nx.spring_layout(hub_net, seed=POS_SEED)  # Seed layout for reproducibility
    nx.draw(hub_net,pos,edge_color=edge_colors['color_list'],node_color=node_colors['color_list'],ax=ax[0],
            with_labels=True, font_weight='bold',font_size=9,node_size=[HUB_NODE_SIZE if v in hub_gene else BASIC_NODE_SIZE for v in hub_net.nodes()])
    fig.colorbar(edge_colors['scale_colormap'], cax=ax[1], orientation='horizontal', label='STRING combined_score')
    fig.colorbar(node_colors['scale_colormap'], cax=ax[2], orientation='horizontal', label='logFC')
    plt.subplots_adjust(hspace = hspace)
    plt.savefig(save)

def _select_grn(tfs,logFC_df,grn,source='source',target='target',pathway_genes=None):
    '''
    Select nodes from input tfs and input pathway

    Parameters
    ----------
    tfs
        Transcription factor names in GRN.
    grn
        Dataframe of gene regulatory network (GRN).
    source
        Column names in grn dataframe with source nodes.
    target
        Column names in grn dataframe with target nodes.
    logFC_df
        Dataframe of deseq2 result.
    pathway_genes
        Genes list of pathway.
    ----------

    Returns
        Dataframe of selected gene regulatory network (GRN)
    '''
    if pathway_genes==None:
        select_grn=grn.loc[[i for i in grn.index if grn.loc[i,source] in tfs],:].reset_index(drop=True)
        select_grn=select_grn.loc[[i for i in select_grn.index if select_grn.loc[i,target] in logFC_df.index],:].reset_index(drop=True)
    else:
        select_grn=grn.loc[[i for i in grn.index if grn.loc[i,source] in tfs],:].reset_index(drop=True)
        pathway_genes = [i for i in pathway_genes if i in logFC_df.index]
        select_grn=select_grn.loc[[i for i in select_grn.index if select_grn.loc[i,target] in pathway_genes],:].reset_index(drop=True)
    return(select_grn)

def tf_net(tfs,grn,logFC_df,name,source='source',target='target',pathway_genes=None,save=TF_outdir+'/TF_net.pdf',
           figsize=(6.8,6),wspace = 0.04,node_cmap='Reds',node_min_weight=None,node_max_weight=None,node_unit=0.1):
    '''
    Plot gene regulatory network (GRN) of input transcription factors. Node color represent logFC, larger Nodes represent input TFs.

    Parameters
    ----------
    tfs
        Transcription factor names in GRN.
    grn
        Gene regulatory network (GRN).
    logFC_df
        dataframe of deseq2 result.
    name
        Column name in logFC_df with gene's logFoldChange.
    source
        Column names in grn dataframe with source nodes.
    target
        Column names in grn dataframe with target nodes.
    pathway_genes
        Genes list of pathway. Optional.
    save
        Path to where to save the plot. Infer the filetype if ending on {`.pdf`, `.png`, `.svg`}.
    ----------

    '''
    select_grn = _select_grn(tfs=tfs,logFC_df=logFC_df,grn=grn,source=source,target=target,pathway_genes=pathway_genes)
    tf_net = nx.from_pandas_edgelist(df=select_grn,source=source,target=target,create_using=nx.DiGraph())
    node_colors = _node_colors(graph=tf_net,logFC_df=logFC_df,name=name,cmap=node_cmap,min_weight=node_min_weight,max_weight=node_max_weight,unit=node_unit)
    fig, ax=plt.subplots(1,2,figsize=figsize,width_ratios=[28,1])
    pos = nx.circular_layout(tf_net)
    nx.draw(tf_net,pos,ax=ax[0],connectionstyle="arc3,rad=0.1",with_labels=True, edge_color='grey',width=0.5,arrowsize=8,
            #edgecolors=['red' if results_df.loc[v,'padj']<0.05 else node_colors['color_list'][list(tf_net.nodes()).index(v)] for v in tf_net.nodes()],
            font_size=9,font_weight='normal',verticalalignment='baseline',
            node_color=node_colors['color_list'],
            node_size=[HUB_NODE_SIZE if v in tfs else BASIC_NODE_SIZE for v in tf_net.nodes()])

    fig.colorbar(node_colors['scale_colormap'], cax=ax[1], orientation='vertical', label='logFC')#horizontal
    plt.subplots_adjust(wspace = wspace)
    plt.savefig(save)

