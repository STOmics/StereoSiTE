#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

"""
@File        :   query_network.py
@Description :   query_network.py
@Author      :   Liuchuandong liuchuandong@genomics.cn
@Date        :   2023/08/01 11:42:48
"""



import os, sys
import pandas as pd
import networkx as nx


####################################################
## 1. query input genes and convert to protein id ##
####################################################
def _query_input(query_term, protein_info, protein_alias):
    '''
    query input genes and convert to protein id
    '''
    query_term = pd.DataFrame(query_term)
    # match preferred_name
    query_term['preferred_name'] = [x if x in list(protein_info['preferred_name']) else None for x in list(query_term[0])]
    # match alias
    query_term['alias'] = [x if x in list(protein_alias['alias']) else None for x in list(query_term[0])]
    # match string_protein_id (which is unique id) according to preferred_name 
    query_term['#string_protein_id'] = [protein_info.loc[list(protein_info['preferred_name']).index(x),'#string_protein_id'] if x!=None else None for x in list(query_term['preferred_name'])]
    # match string_protein_id (which is unique id) according to alias
    rest = query_term[query_term['#string_protein_id'].isna()]
    rest = rest.index
    query_term['#string_protein_id'][rest] = [protein_alias.loc[list(protein_alias['alias']).index(x),'#string_protein_id'] if x!=None else None for x in list(query_term.loc[rest,'alias'])]
    # output
    query_term.rename(columns={query_term.columns[0]:'query_term'})
    return query_term



def _get_all_matched_terms(query_term):
    # remove not match terms
    query_proteins = query_term[query_term['#string_protein_id'].notna()]
    # check duplicated string_protein_id
    dup = query_proteins[query_proteins.duplicated('#string_protein_id',keep=False)]
    if len(dup)>1:
        dup_id = set(dup['#string_protein_id'])
        for i in dup_id:
            tmp=dup[dup['#string_protein_id']==i]
            # keep term which have preferred_name
            if len(tmp['preferred_name'].dropna())==1:
                rm_index = set(tmp.index).difference( set(tmp['preferred_name'].dropna().index) )
                rm_index = list(rm_index)
            # if no preferred_name, keep first in all duplicated
            else:
                rm_index=tmp.index[list( range(1,len(tmp)) )]
            print("remove terms with duplicated protein id: ",flush=True)
            print(query_proteins.loc[rm_index,:],flush=True)
            query_proteins=query_proteins.drop(index=rm_index,axis=0)
    return query_proteins



###################################################
## 2. get STRING network according query term    ##
###################################################
def get_PPInet(query_term, protein_info, protein_alias,full_net,score_cutoff):
    '''
    Get PPI network of query term from STRING database.
    
    Parameters
    ----------
    query_term
        List of input genes.
    protein_info
        Protein information of STRING database.
    protein_alias
        Protein alias of STRING database.
    full_net
        Whole ppi network from STRING database.
    score_cutoff
        STRING combined score of protein-protein interaction, represent confidence level.
    ----------

    Returns
        PPI network of input.
    '''
    query_term = _query_input(query_term, protein_info, protein_alias)
    query_proteins = _get_all_matched_terms(query_term)

    query_net = full_net[full_net['combined_score']>=score_cutoff*1000]
    query_net = query_net[query_net['protein1'].isin(query_proteins['#string_protein_id'])]
    query_net = query_net[query_net['protein2'].isin(query_proteins['#string_protein_id'])]
    query_net = nx.from_pandas_edgelist(query_net,'protein1','protein2',edge_attr = 'combined_score')
    # relabel net
    proteins_id_to_query_dict = query_proteins
    proteins_id_to_query_dict.index = proteins_id_to_query_dict['#string_protein_id']
    proteins_id_to_query_dict = proteins_id_to_query_dict.loc[:,0]
    proteins_id_to_query_dict = proteins_id_to_query_dict.to_dict()
    query_net = nx.relabel_nodes(query_net, proteins_id_to_query_dict, copy=False)
    return query_net


