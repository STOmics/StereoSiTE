#! /usr/bin/env python3
# _*_ coding: utf-8 _*_

#@Author:  LiuXing liuxing2@genomics.cn
#@Date:  2023-07-12 15:29:08
#@Last Modified by: LiuXing
#@Last Modified time: 2023-07-12 15:29:08

import os, sys
import scanpy as sc
import squidpy as sq
import anndata
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import sparse
import math
from tqdm.notebook import tqdm
from multiprocessing import Pool
from optparse import OptionParser
import logging

LOG_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"
logging.basicConfig(level=logging.INFO, format=LOG_FORMAT)

CODE_NUMBER=10000

def _m2h_homologene(adata: anndata):
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

def _mh_translate(adata:anndata):
    if any([x.startswith("Gm") for x in adata.var.index]):
        adata = _m2h_homologene(adata)
        #self.adata.var_names_make_unique()
    return adata

def _hm_translate(adata, LRpairs):
    if "original_gene_symbol" in adata.var.keys():
        #tanslate LRpairs from human genes to mouse genes
        hmdict = adata.var['original_gene_symbol'].to_dict()
        mouse_LRpairs = []
        for LRpair in LRpairs:              
            for gene in LRpair:
                if "_" in gene:
                    mouse_LRpairs.append("_".join([hmdict[z] if z in hmdict.keys() else z for z in gene.split("_")]))
                else:
                    mouse_LRpairs.append(hmdict[gene] if gene in hmdict.keys() else gene)
        return mouse_LRpairs
    else:
        return LRpairs

def _generate_LRpairs(interactionDB:str,
                      adata:anndata,
                      LR_anno:str='annotation'):
    """
    process L-R database.
    """
    file_name = os.path.basename(interactionDB)
    SOURCE = 'source'
    TARGET = 'target'
    if not "MOUSE" in file_name.upper():
        #translate mouse gene symbol to human gene symbol
        adata = _mh_translate(adata)
        logging.info(f"gene homology translation finished.")

    if "CELLCHATDB" in file_name.upper():
        interactions = pd.read_csv(interactionDB)
        interactions[SOURCE] = interactions['ligand']
        interactions[TARGET] = interactions['interaction_name_2'].str.extract('- (.+)', expand=False).map(lambda x: x.strip().strip("(").strip(")").replace("+", "_"))
    else:
        interactions = pd.read_csv(interactionDB)      
        if SOURCE in interactions.columns:
            interactions.pop(SOURCE)
        if TARGET in interactions.columns:
            interactions.pop(TARGET)
        interactions.rename(
                    columns={"genesymbol_intercell_source": SOURCE, "genesymbol_intercell_target": TARGET}, inplace=True
                )
        interactions[SOURCE] = interactions[SOURCE].str.replace("^COMPLEX:", "", regex=True)
        interactions[TARGET] = interactions[TARGET].str.replace("^COMPLEX:", "", regex=True)
    
    LRpairs_dict = dict()
    if LR_anno in interactions.columns:
        LR_types = interactions[LR_anno].unique()
        for LR_type in LR_types:
            LRpairs_dict[LR_type] = interactions[interactions[LR_anno] == LR_type][[SOURCE, TARGET]].drop_duplicates().values
    else:
        LRpairs_dict['All'] = interactions[[SOURCE, TARGET]].drop_duplicates().values
    LRpairs = interactions[[SOURCE, TARGET]].drop_duplicates().values
    LRlist = [] 
    filter_LRpairs = dict()
    for LR_type, LRpairs in LRpairs_dict.items():
        filter_LRpairs[LR_type] = list()
        for LRpair in LRpairs:
            ligand, receptor = LRpair
            ligand_subs = ligand.split("_")
            receptor_subs = receptor.split("_")
            genes = ligand_subs + receptor_subs
            if all(g in adata.var_names for g in genes):
                filter_LRpairs[LR_type].append(LRpair)
                LRlist.extend(genes) 
    adata = adata[:, adata.var_names.isin(LRlist)]
    return filter_LRpairs, adata

def _generate_neighbors_graph(adata:anndata,
                              radius = 400, #int or dict
                              anno:str = "cell2loc_anno"):
        if isinstance(radius, int):
            sq.gr.spatial_neighbors(adata, radius=radius, coord_type="generic", key_added=str(radius))
        elif isinstance(radius, dict):
            for LR_type, r in radius.items():
                sq.gr.spatial_neighbors(adata, radius=r, coord_type='generic', key_added=str(r))
        else:
            raise Exception("radius type must be int or dictionary, but get: {0}".format(type(radius)))
        #sq.gr.nhood_enrichment(adata, cluster_key= anno)
        #np.nan_to_num(adata.uns[f"{anno}_nhood_enrichment"]['zscore'], copy=False)  
        return adata

def preprocess_adata(adata:anndata,
                     interactionDB:str,
                     sample_number:int = 1000000,
                     seed:int = 101,
                     radius = 200, #int or dict
                     scale:float = 0.5,
                     LR_anno:str = "annotation",
                     anno:str = "cell2loc_anno",
                     use_raw:bool = True):

    #adata = anndata.read(adata_file)
    if adata.raw and use_raw:
        adata = adata.raw.to_adata()
    #sampling data with 1000000 cells
    sample_rate  = 1
    if adata.shape[0] > sample_number:
        n_obs = adata.shape[0]
        sample_rate = float(sample_number)/float(n_obs)
        obs_index = adata.obs.index.values.copy()
        np.random.seed(seed) #set the random seed to be 101
        np.random.shuffle(obs_index)
        adata = adata[obs_index[0:sample_number],]
        #sc.pp.subsample(self.adata, n_obs=sample_number)
        logging.info("get subsample data from original data, the subsample data shape is {0}".format(adata.shape))

    #generate LRpairs list from interactions
    LRpairs, adata = _generate_LRpairs(interactionDB, adata, LR_anno=LR_anno)
    logging.info(f"generate LRpairs finished, and get {len(LRpairs)} LRpair")
    #change the radius unit from um to dnb and Check if the LR types given by parameter contians all LR types in the interaction database
    if isinstance(radius, int):
        radius = int(radius/scale)
    elif isinstance(radius, dict):
        for LRtype, r in radius.items():
            radius[LRtype] = int(r/scale)
        if not set(LRpairs.keys()).issubset(set(radius.keys())):
            raise Exception("LR types given by parameter radius: {0} don't contain all LR types in the interaction database: {1}".format(LRpairs.keys(), radius.keys()))
    else:
        raise Exception("The type of radius threshold parameter must be int or dict, but we got: {0}".format(type(radius)))
    adata = _generate_neighbors_graph(adata, radius=radius, anno = anno)
    return adata, LRpairs, sample_rate

def _get_LR_connect_matrix(adata:anndata, 
                           LRpair:list, 
                           connect_matrix:sparse.csr_matrix,
                           complex_process:str = 'mean',
                           distance_coefficient = 0,                      
                           ) -> sparse.coo_matrix:
    ligand, receptor = LRpair    
    if "_" in ligand:
        ligands = ligand.split("_")
        if complex_process.upper() == 'MEAN':
            exp_l = adata[:, ligands].X.mean(axis=1).A1
        elif complex_process.upper() == 'MIN':       
            exp_l = adata[:, ligands].X.min(axis=1).toarray()[:,0]
        else:
            raise Exception("complex process model must be mean or min, but got: {0}".format(complex_process))         
    else:
        exp_l = adata[:, ligand].X.toarray()[:,0]
    if "_" in receptor:
        receptors = receptor.split("_")
        if complex_process.upper() == 'MEAN':
            exp_r = adata[:, receptors].X.mean(axis=1).A1
        elif complex_process.upper() == 'MIN':
            exp_r = adata[:, receptors].X.min(axis=1).toarray()[:,0]
        else:
            raise Exception("complex process model must be mean or min, but got: {0}".format(complex_process))
    else:
        exp_r = adata[:, receptor].X.toarray()[:,0]
    l_rows = np.where(exp_l > 0)[0]
    r_cols = np.where(exp_r > 0)[0]
    sub_connect_matrix = connect_matrix[l_rows,:][:,r_cols].todense()
    dst = np.where(sub_connect_matrix > 0)
    distances = sub_connect_matrix[dst]
    connect_exp_lr = exp_l[l_rows[dst[0]]]*[math.exp(-distance_coefficient*d) for d in distances.A1] + exp_r[r_cols[dst[1]]]
    exp_connect_matrix = sparse.coo_matrix((connect_exp_lr, (l_rows[dst[0]], r_cols[dst[1]])), shape=connect_matrix.shape)
    return exp_connect_matrix

def _get_LR_intensity(exp_connect_matrix:sparse.coo_matrix, 
                      cellTypeIndex:np.array, 
                      cellTypeNumber:int
                      ) -> np.matrix:
    senders = cellTypeIndex[exp_connect_matrix.row]
    receivers = cellTypeIndex[exp_connect_matrix.col]
    interaction_matrix = sparse.csr_matrix((exp_connect_matrix.data, (senders, receivers)), shape=(cellTypeNumber, cellTypeNumber))
    
    return interaction_matrix.todense()
    
def _permutation_test(interaction_matrix:np.matrix,
                      exp_connect_matrix:sparse.coo_matrix,
                      cellTypeIndex:np.array,
                      cellTypeNumber:int,
                      n_perms = 1000,
                      seed = 101
                      ) -> np.array:
    pvalues = np.zeros((cellTypeNumber, cellTypeNumber), dtype=np.int32)
    for i in range(n_perms):
        cellTypeIndex_tmp = cellTypeIndex.copy()
        rs = np.random.RandomState(None if seed is None else i + seed)
        rs.shuffle(cellTypeIndex_tmp)
        interaction_matrix_tmp = _get_LR_intensity(exp_connect_matrix, cellTypeIndex_tmp, cellTypeNumber)
        pvalues += np.where(interaction_matrix_tmp > interaction_matrix, 1, 0)
    pvalues = pvalues/float(n_perms)
    pvalues[interaction_matrix == 0] = None
    return pvalues

def _LRpair_process(adata:anndata,
                    LRpair:np.array,
                    connect_matrix:sparse.coo_matrix,
                    cellTypeIndex:np.array,
                    cellTypeNumber:int,
                    seed:int,
                    n_perms:int,
                    complex_process:str = 'mean',
                    distance_coefficient = 0,
                    )->tuple : #np.array
    exp_connect_matrix = _get_LR_connect_matrix(adata, LRpair, connect_matrix, complex_process = complex_process, distance_coefficient=distance_coefficient)
    interaction_matrix = _get_LR_intensity(exp_connect_matrix, cellTypeIndex, cellTypeNumber)
    pvalues = _permutation_test(interaction_matrix, exp_connect_matrix, cellTypeIndex, cellTypeNumber, seed = seed, n_perms=n_perms)
    return interaction_matrix, pvalues

def _result_combined(result_list:list, #list<np.array>
                     LRpairs:np.array,
                     cell_types:np.array
                     ) -> pd.DataFrame:
    columns = []
    for sender in cell_types:
        for receiver in cell_types:
            columns.append([sender, receiver])
    my_columns = pd.MultiIndex.from_tuples(columns, names=['cluster1', 'cluster2'])
    my_index = pd.MultiIndex.from_tuples(LRpairs, names=["source", "target"])
    values = [np.ravel(x) for x in result_list]
    result_df = pd.DataFrame(np.row_stack(values), index = my_index, columns = my_columns)
    return result_df

def intensities_count(adata:anndata,
                      interactionDB:str,
                      distance_threshold = 200, #int or dict
                      distance_coefficient = 0, #float or dict 
                      scale:float = 0.5,
                      LR_anno:str = "annotation",
                      anno:str = "cell_type",
                      seed:int = 101,
                      n_perms:int = 1000,
                      use_raw:bool = True,
                      jobs:int = 1,
                      complex_process_model:str = 'mean',
                      ) -> dict:
    """
    calculate intensities of interactions between all cell type pairs and ligand receptor pairs.

    Parameters
    ----------
    adata
        anndata
    interactionDB
        file that stores ligand receptor pairs
    distance_threshold
        only cell pairs with the distance shorter than distance_threshold will be connected when construct nearest neighbor graph.
        default=200. The unit is µm
        If the ligand-receptor pairs have been clustered into different types, the distance_threshold can receive a dictionary with 
        LR types and corresponding distance thresholds. e.g: {'Secreted Signaling': 200, 'ECM-Receptor': 200, 'Cell-Cell Contact': 30}
    distance_coefficient
        Consider the distance as one of the factor that influence the interaction intensity using the exponential decay formular: C=C0*e^(-k*d).
        The parameter defines the k value in the formular. Default=0, means distance would not influence the interaction intensity.
        If the ligand-receptor pairs have been clustered into different types, the distance_coefficient can receive a dictionary with
        LR types and corresponding coefficient. e.g: {'Secreted Signaling': 1, 'ECM-Receptor': 0.1, 'Cell-Cell Contact': 0}
    scale
        The distance between adjancent spots, the unit is µm. For Stereo-chip, scale=0.5. default=0.5 
    LR_anno
        The name of the column that contains the LR types annotation information, default=annotation.
    anno
        cell type annotation key, default=cell_type
    seed
        specify seed when randomly shifting cell type labels to do permutation test, which will generate a null distribution
        to calculate p value of each interactions. default=101
    n_perms
        specify permutation number. default=10000
    use_raw
        bool value, which ditermine whether use the raw data in anndata. default=True
    jobs
        when jobs > 1, the program will call multi-process to analyze data. jobs=1 is fast enough for most task. default=1
    complex_process_model
        determine how to deal with the complexed ligand and receptor which contain multiple subunits. There are two options: mean, min.
        mean: calculate the mean expression of all subunits to represent the complex
        min: pick the minimal expression of all subunits to represent the complex

    Returns
    ----------
    dictionary contains intensities and p values result
        {'intensities': intensities in DataFrame,
         'pvalues': pvalues in DataFrame
        }
    """

    adata, LRpairs, sample_rate = preprocess_adata(adata, interactionDB, radius=distance_threshold, scale=scale, seed = seed, LR_anno=LR_anno, anno=anno, use_raw=use_raw)
    #connect_matrix = adata.obsp['spatial_connectivities']
    cell_types = adata.obs[anno].unique()
    cell_type_dict = dict(zip(cell_types, range(0, len(cell_types))))
    cellTypeIndex = adata.obs[anno].map(cell_type_dict).astype(int).values
    cellTypeNumber = len(cell_types)

    logging.info("interaction intensity count begin.")
    results = []
    intensities_list = []
    pvalues_list = []
    LRpairs_lists = []
    if (jobs == 1):
        for LRtype, LRpair_list in LRpairs.items():
            logging.info("compute the interaction intensity of LRpairs {0}".format(LRtype))
            if isinstance(distance_threshold, (int, float)):
                connect_matrix = adata.obsp[f"{int(distance_threshold/scale)}_distances"]
            elif isinstance(distance_threshold, dict):
                connect_matrix = adata.obsp[f"{int(distance_threshold[LRtype])}_distances"]
            else:
                raise Exception("the type of distance_threshold must be int, float or dict, but get:{0}".format(type(distance_threshold)))
            if isinstance(distance_coefficient, (int, float)):
                k = distance_coefficient
            elif isinstance(distance_coefficient, dict):
                k = distance_coefficient[LRtype]
            else:
                raise Exception("the type of distance_threshold must be int, float or dict, but get:{0}".format(type(distance_coefficient)))
            for LRpair in LRpair_list:
                LRpairs_lists.append(LRpair)
                results.append(_LRpair_process(adata, LRpair, connect_matrix, cellTypeIndex,cellTypeNumber, seed, n_perms, complex_process=complex_process_model, distance_coefficient=k))
    else:
        pool = Pool(jobs)
        for LRtype, LRpair_list in LRpairs.items():
            logging.info("compute the interaction intensity of LRpairs {0}".format(LRtype))
            if isinstance(distance_threshold, int):
                connect_matrix = adata.obsp[f"{int(distance_threshold/scale)}_distances"]
            elif isinstance(distance_threshold, dict):
                connect_matrix = adata.obsp[f"{int(distance_threshold[LRtype])}_distances"]
            else:
                raise Exception("the type of distance_threshold must be int or dict, but get:{0}".format(type(distance_threshold)))
            if isinstance(distance_coefficient, (int, float)):
                k = distance_coefficient
            elif isinstance(distance_coefficient, dict):
                k = distance_coefficient[LRtype]
            else:
                raise Exception("the type of distance_threshold must be int, float or dict, but get:{0}".format(type(distance_coefficient)))
            LRpairs_lists.extend(LRpair_list)
            for LRpair in tqdm(LRpair_list):
                LRpairs_lists.append(LRpair)
                results.append(pool.apply_async(_LRpair_process, (adata, LRpair, connect_matrix, cellTypeIndex, cellTypeNumber, seed, n_perms, complex_process_model, k,)))
        pool.close()
        pool.join()
    logging.info("interaction intensity count finished.")
    LRpairs_lists = _hm_translate(adata, LRpairs_lists)
    logging.info("begin to combine results")
    for result in results:
        if (jobs == 1):
            intensities, pvalues = result
        else:
            intensities, pvalues = result.get()
        intensities_list.append(intensities)
        pvalues_list.append(pvalues)
    intensity_df = _result_combined(intensities_list, LRpairs_lists, cell_types)
    pvalues_df = _result_combined(pvalues_list, LRpairs_lists, cell_types)
    intensity_df = intensity_df/sample_rate
    logging.info(f"result combining finished.")
    
    plot_data = {'intensities': intensity_df,
                'pvalues': pvalues_df}
    return plot_data

def intensities_write(plot_data, out_dir):

    if not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    intensity_df = plot_data['intensities']
    pvalues_df = plot_data['pvalues']
    intensity_df_pickle = f"{out_dir}/intensities.pkl"     
    pvalues_df_pickle = f"{out_dir}/pvalues.pkl"      
    intensity_df.to_pickle(intensity_df_pickle)    
    pvalues_df.to_pickle(pvalues_df_pickle) 
    logging.info("finished to write result to pickle file.")

    intensity_df_csv = f"{out_dir}/intensities.csv"
    pvalues_df_csv = f"{out_dir}/pvalues.csv"
    intensity_df.columns = ["|".join(x) for x in intensity_df.columns]
    pvalues_df.columns = ["|".join(x) for x in pvalues_df.columns]
    intensity_df.to_csv(intensity_df_csv)
    pvalues_df.to_csv(pvalues_df_csv)
    logging.info("finished to write result to csv file.")

def intensities_with_radius(adata: anndata, 
                            plot_pairs_file: str,
                            radius_list: list = list(range(10, 400, 10)),
                            anno: str = "cell_type",
                            copy=False):
    """
    calculate the intensities with different radius threshold

    Parameters
    ----------
    adata
        anndata
    plot_pairs_file
        path of csv file which stores the interactions that will be processed.
        Contains 5 columns: sender, receiver, ligand, receptor. Columns was seperated by table
    radius_list
        list contains radius thresholds that will be set to calculate intensity
    anno
        cell type annotation key
    copy
        If true, result will be return in DataFrame format. If false, result will be stored in anndata.uns['intensities_with_radius'].
    """
    radius_list = [x*2 for x in radius_list]
    cell_types = adata.obs[anno].unique()
    cell_type_dict = dict(zip(cell_types, range(0, len(cell_types))))
    cellTypeIndex = adata.obs[anno].map(cell_type_dict).astype(int).values
    cellTypeNumber = len(cell_types)
    
    plot_pairs = pd.read_csv(plot_pairs_file, sep="\t")
    genes_list = plot_pairs[['ligand', 'receptor']].drop_duplicates().values.tolist()
    results = {}
    for radius in tqdm(radius_list):
        connectivities_key = f"{radius}_connectivities"
        if not connectivities_key in adata.obsp.keys():
            connect_matrix = sq.gr.spatial_neighbors(adata, radius=radius, coord_type="generic", copy=True)
        else:
            connect_matrix = adata.obsp[connectivities_key]
        intensities_list = []
        for genes in genes_list:          
            exp_connect_matrix = _get_LR_connect_matrix(adata, genes, connect_matrix)
            interaction_matrix = _get_LR_intensity(exp_connect_matrix, cellTypeIndex, cellTypeNumber)
            intensities_list.append(interaction_matrix)
        intensity_df = _result_combined(intensities_list, genes_list, cell_types)
        results[radius] = intensity_df
    plot_df = pd.DataFrame()
    plot_df['radius'] = [x/2 for x in radius_list]
    for index, pair in plot_pairs.iterrows():
        intensities = []
        for radius in radius_list:
            intensities.append(results[radius].loc[tuple(pair[['ligand', 'receptor']].values)][tuple(pair[['sender', 'receiver']].values)])
        column = f"{pair['sender']}|{pair['receiver']} ({pair['ligand']}|{pair['receptor']})"
        plot_df[column] = intensities
    if copy:
        return plot_df
    else:
        adata.uns['intensities_with_radius'] = plot_df

def main():
    """
    This program can calculate the spatial cell interaction intensity
    %prog [options]
    """
    parser = OptionParser(main.__doc__)
    parser.add_option("-i", "--adata", action = "store", type = "str", dest = "adata", help = "input stereo-seq data in h5ad format.")
    parser.add_option("-o", "--out_dir", action = "store", type = "str", dest = "out_dir", help = "output directory path.")
    parser.add_option("-d", "--interactionDB", action = "store", type = "str", dest = "interactions", help = "interaction database file path.")
    parser.add_option("--distances", action = "store", type = "str", dest = "distances", default = "all=200", help = "if the distance between adjacent cells is less than the threshold, they will be connected.\
                       The unit is µm. default: all=200 (µm), or you can give different distance threshold for individual LR types, such as: 'Secreted Signaling=200,ECM-Receptor=200,Cell-Cell Contact=30'")
    parser.add_option("--distances_coefficient", action="store", type = "str", dest="distances_coefficient", default = "all=0", help = "If the distances_coefficient was not equal to 0, when calculate the \
                      interaction intensity between sender and receiver cells, the distances would be considered using the exponential decay formular: C=C0*e^(-k*d). The parameter defines the k value in the formular.\
                      default: all=0, or you can give different distance coefficient for individual LR types, such as: 'Secreted Signaling=1,ECM-Receptor=0.1,Cell-Cell Contact=0'."
                      )
    parser.add_option("-t", "--threads", action = "store", type = "int", dest = "threads", default = 1, help = "number of processing that will be used to run this program. default=1")
    parser.add_option("--n_perms", action = "store", type = "int", dest = "n_perms", default = 1000, help = "number of permutation test times. default=1000")
    parser.add_option("--seed", action = "store", type = "int", dest = "seed", default = 101, help = "seed for randomly shuffle cell type label. default=101")
    parser.add_option("--use_raw", action = "store", type = "str", dest = "use_raw", default = "True", help = "if use raw data in the anndata. default=True")
    parser.add_option("--cell_type", action = "store", type = "str", dest = "cell_type", default = "cell2loc_anno", help = "cell type annotation key. default=cell2loc_anno")
    parser.add_option("--complex_process", action = "store", type = "str", dest = "complex_process", default = "mean", help = "choose the procession method for complex, mean or min. default=mean")
    opts, args = parser.parse_args()

    if (opts.adata == None or opts.out_dir == None or opts.interactions == None):
        sys.exit(not parser.print_help())
    distances_str = opts.distances
    distances = dict([[x.strip().split("=")[0], float(x.strip().split("=")[1])] for x in distances_str.strip().split(",")])
    distances_coefficient_str = opts.distances_coefficient
    distances_coefficient = dict([[x.strip().split("=")[0], float(x.strip().split("=")[1])] for x in distances_coefficient_str.strip().split(",")])
    if 'all' in distances.keys():
        distances = distances['all']
    if 'all' in distances_coefficient.keys():
        distances_coefficient = distances_coefficient['all']
    adata = anndata.read(opts.adata)
    plot_data = intensities_count(adata, 
                      opts.interactions, 
                      distance_threshold=distances,
                      distance_coefficient=distances_coefficient,
                      anno=opts.cell_type,
                      n_perms=opts.n_perms, 
                      seed=opts.seed,
                      jobs = opts.threads,
                      use_raw=bool(opts.use_raw),
                      complex_process_model = opts.complex_process)
    intensities_write(plot_data, opts.out_dir)

if __name__ == "__main__":
    main()
        