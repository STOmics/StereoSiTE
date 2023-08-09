#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

"""
@File        :   ppi_analysis.py
@Description :   ppi_analysis.py
@Author      :   Liuchuandong liuchuandong@genomics.cn
@Date        :   2023/08/01 14:35:03
"""

from gc import DEBUG_COLLECTABLE
import os, sys
from optparse import OptionParser
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import subprocess
import logging

PPI_outdir = 'PPI'
os.makedirs(PPI_outdir, exist_ok=True)
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def mcl(query_net,inflation: float = 4.0):
    '''
    Run Markov CLustering in query net

    Parameters
    ----------
    query_net
        the ppi net queried from STRING database according to user input
    inflation
        the Markov CLuster algorithm's parparmeter, vary this parameter to obtain clusterings at different levels of granularity. A good set of starting values is 1.4, 2, 4, and 6.
    ----------
    '''
    logging.critical("Make sure the Markov CLuster algorithm (MCL,https://github.com/micans/mcl) be installed properly.")
    
    nx.to_pandas_edgelist(query_net).to_csv(PPI_outdir+"/query_net_edgelist.txt",sep='\t',header=0,index=0)
    mcl_result_filepath = PPI_outdir+'/out.query_net_edgelist.txt.I'+str( int(inflation*10) )
    
    subprocess.run(['mcl', PPI_outdir+'/query_net_edgelist.txt','--abc','-I', str(inflation),'-o',mcl_result_filepath], capture_output=True)
    logging.info('mcl result be stored in: '+mcl_result_filepath )
    mcl_result = pd.read_csv(mcl_result_filepath, sep='\t', header=None, index_col=False)

    logging.info('Get {} clusters by markov clustering'.format( len(mcl_result) ))
    logging.info('The lagrest cluster with {} nodes'.format( len(mcl_result.loc[0,].dropna()) ))

    return(mcl_result)



def get_cluster_net(query_net,mcl_result,mcl_id=0):
    '''
    Subset specific net from qurey net according to specific mcl_cluster_id, default is the largest one 
    '''
    cluster_net = nx.subgraph(query_net,[i for i in list(mcl_result.loc[mcl_id,]) if str(i)!='nan'])
    return(cluster_net)



def get_MCC_hub_genes(cluster_net):
    '''
    Find hub genes from specific MCL cluster
    '''
    # 1. get genes by sorting MCC score
    node_MCC  = nx.node_clique_number(cluster_net)
    rank_values = set(node_MCC.values())
    rank_values = sorted(rank_values,reverse=True)
    max_values = rank_values[0:10]
    print('Top10 MCC scores are:',max_values)
    max_MCC_genes = []
    # get genes which with top MCC score, by default output 10 genes, if more than 10 genes with max score, still output
    for i in max_values:
        for k,v in node_MCC.items():
            if v ==i: max_MCC_genes.append(k)
        if len(max_MCC_genes)>10:
            break
    #  2. get genes by sorting nodes degree
    top_degree_genes=[]
    rank_degree = sorted(cluster_net.degree, key=lambda x: x[1], reverse=True)
    degree = list(set([i[1] for i in rank_degree]))
    degree.sort(reverse=True)
    if len(degree)<10:
        degree_cutoff=degree[len(degree)-1]
    else:
        degree_cutoff=degree[9]
    for i in rank_degree:
        if i[1]>=degree_cutoff:
            top_degree_genes.append(i[0])
    print('Top10 degree are:',degree[0:10])
    # 3. get overlapped genes according to MCC and degree
    hub_genes = [i for i in max_MCC_genes if i in top_degree_genes]
    hub_genes_df = pd.DataFrame({'hub_genes':hub_genes})
    rank_degree = dict(rank_degree)
    hub_genes_df['degree'] = [rank_degree[i] for i in hub_genes]
    hub_genes_df['MCC_score'] = [node_MCC[i] for i in hub_genes]
    hub_genes_df = hub_genes_df.sort_values(by='degree',ascending=False,ignore_index=True)
    return(hub_genes_df)



def get_hub_net(hub_genes,cluster_net,cutoff: float = 0.8):
    '''
    Subset hub net of hub gene from specific MCL cluster

    Parpameters
    hub_genes
        The hub genes from specific MCL cluster by using get_hub_gene function.
    cluster_net
        The specific MCL cluster.
    cutoff
        the edge's STRING combined-score cutoff.
    '''
    # 1. remove edges which score<cutoff
    cluster_net=cluster_net.copy()
    for i in list(cluster_net.edges()):
        if cluster_net.edges[i]['combined_score'] < cutoff*1000:
            cluster_net.remove_edge(i[0],i[1])
    # 2. get hub genes and their neighbors
    hub_genes_set=set(hub_genes)
    for i in hub_genes:
        neigb = [n for n in cluster_net.neighbors(i)]
        hub_genes_set = hub_genes_set.union(set(neigb))
    hub_net = nx.subgraph(cluster_net,list(hub_genes_set))
    return(hub_net)

