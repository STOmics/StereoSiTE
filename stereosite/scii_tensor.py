import numpy as np
import pandas as pd
import squidpy as sq
import os, sys
from itertools import product
from scipy import sparse
import time
import tensorly as tl
from tensorly.decomposition import non_negative_tucker, non_negative_tucker_hals
from tqdm import tqdm
import logging as logg
import torch
import random

from stereosite.scii import _generate_LRpairs, _get_LR_connect_matrix, _get_LR_intensity, _permutation_test,  _result_combined, _mh_translate


format="%(asctime)s-%(levelname)s-%(message)s"
logg.basicConfig(format=format, level=logg.INFO)

def get_keys(d, value):
    return [k for k, v in d.items() if v == value]

def _sorted(s, num):
    tmp = s.sort_values(ascending=False)[:num].index
    tmp.index = range(num)
    return tmp

def find_max_column_indices(array):
    max_column_indices = np.argmax(array, axis=1)
    return max_column_indices

class InteractionTensor():
    """
    Class for constructing and analyzing Spatial Cellular Interaction Intensity (SCII) matrices
    based on spatial transcriptomics data.

    The InteractionTensor class provides methods to build SCII matrices, process them,
    evaluate tensor ranks, and perform tensor factorization for further analysis.

    Attributes:
    - adata: anndata
        The spatial transcriptomics data.
    - interactionDB: str
        Database of ligand-receptor pairs.
    """
    
    def __init__(self, adata, interactionDB: str=None):
        # Save variables for this class
        self.adata = adata
        self.interactionDB = interactionDB

def build_SCII(interactiontensor: InteractionTensor=None, radius = 200, #int or dict 
               scale:float=0.5,
               coord_type: str="generic",
               window_size: int=400, anno_col: str='cell2loc_anno', use_raw: bool=True, interactionDB: str=None) -> list:
    """
    Construct a Spatial Cellular Interaction Intensity (SCII) matrix based on spatial transcriptomics data.

    This function calculates the interactions between different cell types within specified spatial windows
    and generates a CCI matrix containing interaction strengths.

    Parameters
    ----------
    interactiontensor : InteractionTensor
        The InteractionTensor of the spatial transcriptomics data.
    radius : int
        The radius for spatial neighbor search, by default 200 µm. 
        If the ligand-receptor pairs have been clustered into different types, the distance_threshold can receive a dictionary with 
        LR types and corresponding distance thresholds. e.g: {'Secreted Signaling': 200, 'ECM-Receptor': 200, 'Cell-Cell Contact': 30}
    scale : float
        The distance between adjancent spots, the unit is µm. For Stereo-chip, scale=0.5. default=0.5
    coord_type : str
        The coordinate type for spatial neighbor calculation, by default "generic".
    window_size : int
        The size of the spatial windows, by default 400µm.
    anno_col : str
        Column name containing cell annotations, by default 'cell2loc_anno'.
    use_raw : bool
        whether use raw counts to build CCI matrix, by default True.
    interactionDB : str
        Database of ligand-receptor pairs.

    Returns
    -------
    list of pandas.DataFrame
        A list of DataFrames containing interaction strengths between cell types in each window.
    """
    interactionDB = interactiontensor.interactionDB
    adata = interactiontensor.adata

    interactiontensor.radius = radius
    interactiontensor.scale = scale
    interactiontensor.coord_type = coord_type
    interactiontensor.window_size = window_size
    interactiontensor.anno_col = anno_col
    interactiontensor.use_raw = use_raw

    time_start=time.time()

    if use_raw and adata.raw != None:
        adata = adata.raw.to_adata()
    #else:
    #    sc.pp.normalize_total(adata, inplace=True, target_sum=1e4)
    #    sc.pp.log1p(adata)

    # filter LR genes
    logg.info("Filter LR genes")
    filter_LRpairs_dict, adata_sub = _generate_LRpairs(interactionDB=interactionDB, adata=adata)

    # nearest neighbor graph
    logg.info("Create spatial neighbor graph")
    if isinstance(radius, int):
        radius = radius/scale
        sq.gr.spatial_neighbors(adata_sub, radius=radius, coord_type=coord_type, key_added=str(radius))
    elif isinstance(radius, dict):
        for LRtype, radii in radius.items():
            radii = radii/scale
            radius[LRtype] = radii
            sq.gr.spatial_neighbors(adata_sub, radius=radii, coord_type=coord_type, key_added=str(radii))
    else:
        raise Exception("radius type must be int or dict, but we got {0}".format(type(radius)))
 
    exp_connect_list = []
    filter_LRpairs = []
    for LRtype, LRpair_list in filter_LRpairs_dict.items():
        if isinstance(radius, int):
            connect_matrix = adata_sub.obsp[f"{radius}_connectivities"]
        elif isinstance(radius, dict):
            if not LRtype in radius.keys():
                raise Exception(f"The radius threshold of LRtype {LRtype} in the interactionDB was not specified by the radius parameter: {radius}")
            connect_matrix = adata_sub.obsp[f"{radius[LRtype]}_connectivities"]
        for genes in LRpair_list:
            filter_LRpairs.append(genes)
            exp_connect_matrix = _get_LR_connect_matrix(adata_sub, genes, connect_matrix)
            exp_connect_list.append(exp_connect_matrix)
    interactiontensor.filter_LRpair = filter_LRpairs

    interactiontensor.exp_connect_list = exp_connect_list

    window_size = window_size/scale
    adata_sub.obs[f'x_{window_size}'] = (adata_sub.obsm['spatial'][:, 0]//window_size)*window_size+(window_size/2)
    adata_sub.obs[f'y_{window_size}'] = (adata_sub.obsm['spatial'][:, 1]//window_size)*window_size+(window_size/2)
    adata_sub.obs[f'bin{window_size}'] = adata_sub.obs[f'x_{window_size}'].astype('str') + '_' + adata_sub.obs[f'y_{window_size}'].astype("str")

    adata_sub.obs['cellname'] = adata_sub.obs.index.tolist()
    df = adata_sub.obs[['cellname', f'bin{window_size}']].reset_index()
    df['index'] = df.index.tolist()
    bin_group = df[['index', f'bin{window_size}']].groupby(f'bin{window_size}')
    indices = [indices for t, indices in bin_group.groups.items()]

    interactiontensor.indices = indices

    logg.info("Start build SCII matrix")
    results = []
    # i=0
    celltype = adata_sub.obs[anno_col].values.unique().tolist()
    combinations = list(product(celltype, repeat=2))
    for i in tqdm(range(len(indices))):   
        store_mt = pd.DataFrame(np.zeros((len(combinations), len(filter_LRpairs))), index = ['_'.join(comb) for comb in combinations], columns=["-".join(x) for x in filter_LRpairs]).astype('float32')

        adata_sub_sub = adata_sub[indices[i], :]

        cell_types = adata_sub_sub.obs[anno_col].unique()
        cell_type_dict = dict(zip(cell_types, range(0, len(cell_types))))
        cellTypeIndex = adata_sub_sub.obs[anno_col].map(cell_type_dict).astype(int).values

        cellTypeNumber = len(adata_sub_sub.obs[anno_col].unique())

        for j in range(len(filter_LRpairs)):
            # j=0
            lr = "-".join(filter_LRpairs[j])

            data = exp_connect_list[j].tocsr()[indices[i], :][:, indices[i]].tocoo()

            senders = cellTypeIndex[data.row]
            receivers = cellTypeIndex[data.col]
            interaction_matrix = sparse.csr_matrix((data.data, (senders, receivers)), shape=(cellTypeNumber, cellTypeNumber))

            nonzero_row, nonzero_col = interaction_matrix.nonzero()
            for r in range(len(nonzero_row)):
                for c in range(len(nonzero_col)):
                    cellpair = get_keys(cell_type_dict, nonzero_row[r])[0] + "_" + get_keys(cell_type_dict, nonzero_col[c])[0]

                    store_mt.loc[cellpair][lr] = interaction_matrix[nonzero_row[r], nonzero_col[c]]
        results.append(store_mt)

    interactiontensor.lr_mt_list = results
    time_end=time.time()
    logg.info(f"Finish build CCI matrix - time cost {(((time_end - time_start) / 60) / 60)} h")


def process_SCII(interactiontensor: InteractionTensor=None, 
            zero_remove: bool=True,
            log_data: bool=True) -> np.ndarray:

    """
    Process the Spatial Cellular Interaction Intensity (SCII) matrix list.

    Parameters
    ----------
    interactiontensor : InteractionTensor
        The InteractionTensor of the spatial transcriptomics data.
    bin_zero_remove : bool
        Flag indicating whether to remove spatial bins with zero intensity, by default True.
    log_data : bool
        Flag indicating whether to log-transform the data, by default True.

    Returns
    -------
    np.ndarray
        Processed 3D Celltype-Celltype Interaction (CCI) matrix.
    """

    time_start=time.time()
    lr_mt_list = interactiontensor.lr_mt_list
    final_mt = np.dstack(lr_mt_list)

    if zero_remove == True:
        window_zero_indices = np.all(final_mt == 0, axis=(0,1))
        lrpair_zero_indices = np.all(final_mt == 0, axis=(0,2))
        cellpair_zero_indices = np.all(final_mt == 0, axis=(1,2))
        
        logg.info(f"{sum(window_zero_indices)} window, {sum(lrpair_zero_indices)} lrpair, {sum(cellpair_zero_indices)} cellpair have zero intensity")
        
        interactiontensor.zero_indices = [cellpair_zero_indices, lrpair_zero_indices, window_zero_indices]
        
        final_mt = final_mt[:, :, ~window_zero_indices]
        final_mt = final_mt[:, ~lrpair_zero_indices, :]
        final_mt = final_mt[~cellpair_zero_indices, :, :]
        
        cellpair = lr_mt_list[0].index[~cellpair_zero_indices]
        interactiontensor.cellpair = cellpair
    
        lrpair = lr_mt_list[0].columns[~lrpair_zero_indices]
        interactiontensor.lrpair = lrpair
        
        lr_mt_list_filter = []
        for i in range(len(lr_mt_list)):
            lr_mt_list_filter.append(lr_mt_list[i].loc[cellpair][lrpair])
        lr_mt_list_filter = [k for k,v in zip(lr_mt_list_filter, ~window_zero_indices) if v]
        interactiontensor.lr_mt_list_filter = lr_mt_list_filter
        
        indices = interactiontensor.indices

        indices_filter = [k for k,v in zip(indices, (~window_zero_indices).tolist()) if v]
        interactiontensor.indices_filter = indices_filter

    else:
        interactiontensor.cellpair = lr_mt_list[0].index
        interactiontensor.lrpair = lr_mt_list[0].columns
        interactiontensor.zero_indices = []
        interactiontensor.lr_mt_list_filter = lr_mt_list
        interactiontensor.indices_filter = indices

    if log_data:
        final_mt = np.log1p(final_mt)

    interactiontensor.cci_matrix = final_mt

    time_end=time.time()
    logg.info(f"Finish processing LR matrix - time cost {(((time_end - time_start) / 60) / 60)} h")

def evaluate_ranks(interactiontensor: InteractionTensor=None, 
                   num_TME_modules: int=20, num_cellpair_modules: int=20, 
                   device: str="cuda:0", init: str='svd', 
                   n_iter_max: int=10, 
                   method: str='tucker', 
                   use_gpu: bool=True,
                   backend:str='pytorch'
                   ) -> np.ndarray:
    """
    Evaluate reconstruction error of tensor decomposition of Cell-Cell Interaction (CCI) tensor matrix.

    Parameters
    ----------
    interactiontensor : InteractionTensor
        The InteractionTensor of the spatial transcriptomics data.
    num_TME_modules : int, optional
        Number of modules for the TME (Tumor Microenvironment) to test, by default 20.
    num_cellpair_modules : int, optional
        Number of modules for cell-cell pairs and ligand-receptor pairs to test, by default 20.
    device : str, optional
        Device for computation ('cpu' or 'cuda:X'), by default 'cuda:0'.
    init : str, optional
        Initialization method for tensor decomposition, by default 'svd'.
    n_iter_max : int, optional
        Maximum number of iterations for tensor decomposition, by default 10.
    method : str, optional
        Method for tensor decomposition ('tucker' or 'hals'), by default 'tucker'.
    use_gpu : bool, optional
        Whether to use GPU for computation, by default True.

    Returns
    -------
    np.ndarray
        A matrix representing evaluated reconstruction errors.
    """

    time_start=time.time()
    tl.set_backend(backend)
    dat = interactiontensor.cci_matrix
    torch.cuda.empty_cache()
    # pal = sns.color_palette('bright',10)
    # palg = sns.color_palette('Greys',10)
    num_TME_modules = num_TME_modules + 1
    mat1 = np.zeros((num_TME_modules,num_cellpair_modules))
    # mat2 = np.zeros((num_TME_modules,num_cellpair_modules))
    if use_gpu:
        tensor = tl.tensor(dat, device=device)
    else:
        tensor = tl.tensor(dat)
        tensor = tensor.to("cpu")

    for i in tqdm(range(2,num_cellpair_modules)):
        for j in range(1,num_TME_modules):
            # we use NNTD as described in the paper

            if method == 'tucker':
                facs_overall = non_negative_tucker(tensor,rank=[i,i, j],random_state = 2337, init=init, n_iter_max=n_iter_max)
                # facs_overall = [factor.cpu() for factor in facs_overall]
                mat1[j,i] = np.mean((dat- tl.to_numpy(tl.tucker_to_tensor((facs_overall[0],facs_overall[1]))))**2)
                # mat2[j,i] = np.linalg.norm(dat - tl.to_numpy(tl.tucker_to_tensor((facs_overall[0],facs_overall[1]))) ) / np.linalg.norm(dat)
            elif method == 'hals':
                facs_overall = non_negative_tucker_hals(tensor,rank=[i,i, j],random_state = 2337, init=init, n_iter_max=n_iter_max)
                mat1[j,i] = np.mean((dat- tl.to_numpy(tl.tucker_to_tensor((facs_overall[0],facs_overall[1]))))**2)
    interactiontensor.reconstruction_errors = mat1
    time_end=time.time()
    logg.info(f"Finish eval SCII tensor rank - time cost {(((time_end - time_start) / 60) / 60)} h")
    return mat1

def SCII_Tensor(interactiontensor: InteractionTensor=None, rank: list=[8,8,8], random_state: int=32, init: str="svd", n_iter_max: int=100, normalize_factors:bool=True,
                backend: str="pytorch",
                device='cuda:0'):
    """
    Perform tensor factorization and analysis on a Cell-Cell Interaction (CCI) matrix.

    This function reads CCI matrix data, performs non-negative tensor factorization, generates visualizations,
    and saves the results for further analysis.

    Parameters
    ----------
    interactiontensor : InteractionTensor
        The InteractionTensor of the spatial transcriptomics data.
    rank : list, optional
        Rank for the non-negative tensor factorization, by default [8, 8, 8].
    random_state : int, optional
        Random state for reproducibility, by default 32.
    init : str, optional
        Initialization method for tensor factorization, by default 'svd'.
    n_iter_max : int, optional
        Maximum number of iterations for tensor factorization, by default 100.
    backend : str, optional
        Backend for tensor operations, by default 'pytorch'.
    device : str, optional
        Device for computation ('cpu' or 'cuda:X'), by default 'cuda:0'.
    """

    time_start=time.time()

    torch.cuda.empty_cache()
    tl.set_backend(backend)

    tensor = tl.tensor(interactiontensor.cci_matrix, device=device)
    core, factors = non_negative_tucker(tensor, rank=rank, random_state = random_state, init=init, n_iter_max=n_iter_max, normalize_factors=normalize_factors)

    factors = [x.data.cpu().numpy() for x in factors]
    core = core.data.cpu().numpy()

    interactiontensor.factors = factors
    interactiontensor.core = core

    cellpair = interactiontensor.cellpair
    lrpair = interactiontensor.lrpair

    tme_cluster = find_max_column_indices(factors[2])
    interactiontensor.tme_cluster = tme_cluster

    indices_filter = interactiontensor.indices_filter

    adata = interactiontensor.adata
    adata.obs['TME_module'] = None
    for k,v in enumerate(tme_cluster):
        adata.obs['TME_module'].iloc[indices_filter[k]] = str(v)

    adata.obs.TME_module = adata.obs.TME_module.astype('category')
    adata.obs.TME_module = adata.obs.TME_module.cat.set_categories(np.arange(adata.obs.TME_module.cat.categories.astype('int').max()+1).astype('str'))
    interactiontensor.adata = adata

    lr_df = pd.DataFrame(factors[1], index = lrpair)
    interactiontensor.lr_factor = lr_df

    cc_df = pd.DataFrame(factors[0], index = cellpair)
    interactiontensor.cc_factor = cc_df

    time_end=time.time()
    logg.info(f"Finish SCII tensor - time cost {(((time_end - time_start) / 60) / 60)} h")


def SCII_Tensor_multiple(interactiontensor: InteractionTensor=None, rank: list=[8,8,8], random_state: int=32, init: str="svd", n_iter_max: int=100, backend: str="pytorch",
                        device: str='cuda:0'):
    """
    Perform tensor factorization and analysis on a Cell-Cell Interaction (CCI) matrix.

    This function reads CCI matrix data, performs non-negative tensor factorization, generates visualizations,
    and saves the results for further analysis.

    Parameters
    ----------
    interactiontensor : InteractionTensor
        The InteractionTensor of the spatial transcriptomics data.
    rank : list, optional
        Rank for the non-negative tensor factorization, by default [8, 8, 8].
    random_state : int, optional
        Random state for reproducibility, by default 32.
    init : str, optional
        Initialization method for tensor factorization, by default 'svd'.
    n_iter_max : int, optional
        Maximum number of iterations for tensor factorization, by default 100.
    backend : str, optional
        Backend for tensor operations, by default 'pytorch'.
    top_n_cc : int, optional
        Number of top cell-cell pairs to consider for visualization, by default 3.
    top_n_lr : int, optional
        Number of top LR pairs to consider for visualization, by default 3.
    figsize : tuple, optional
        Size of the figure for heatmaps, by default (8,15).
    device : str, optional
        Device for computation ('cpu' or 'cuda:X'), by default 'cuda:0'.
    """

    time_start=time.time()

    torch.cuda.empty_cache()
    tl.set_backend(backend)
    tensor = tl.tensor(interactiontensor.cci_matrix, device=device)

    core, factors = non_negative_tucker(tensor, rank=rank, random_state = random_state, init=init, n_iter_max=n_iter_max)

    factors = [x.data.cpu().numpy() for x in factors]
    core = core.data.cpu().numpy()

    interactiontensor.factors = factors
    interactiontensor.core = core

    cellpair = interactiontensor.cellpair
    lrpair = interactiontensor.lrpair

    num_TME_modules=rank[2]

    tme_cluster = find_max_column_indices(factors[2])
    interactiontensor.tme_cluster = tme_cluster
 
    patient_id = list(interactiontensor.module_lr_mt.keys())
    module_lr_mt = interactiontensor.module_lr_mt
    
    individual_patient = [[k]*v.shape[2] for k,v in zip(patient_id, module_lr_mt.values())]
    individual_patient = [item for sublist in individual_patient for item in sublist]
    interactiontensor.patient_map = individual_patient 

    individual_module = [np.arange(v.shape[2]) for v in module_lr_mt.values()]
    individual_module = [item for sublist in individual_module for item in sublist]
    interactiontensor.module_map = individual_module

    adata_list = interactiontensor.adata

    adata_list[0].obs['tmp'] = [str(random.randint(0, num_TME_modules-1)) for _ in range(adata_list[0].n_obs)]

    ids = patient_id
    for i in range(len(ids)):

        tme_cluster_idx = [v for k,v in zip(individual_patient, tme_cluster) if k == ids[i]]

        tme_module_idx = [v for k,v in zip(individual_patient, individual_module) if k == ids[i]]

        adata_list[i].obs['TME_meta_module'] = None

        adata_list[i].obs['TME_meta_module'] = adata_list[i].obs['TME_module'].map(dict(zip([str(x) for x in tme_module_idx], [str(x) for x in tme_cluster_idx])))

        adata_list[i].obs['TME_meta_module'] =  adata_list[i].obs['TME_meta_module'].astype('category')

    interactiontensor.adata_list_TME = adata_list

    lr_df = pd.DataFrame(factors[1], index = lrpair)
    interactiontensor.lr_factor = lr_df

    cc_df = pd.DataFrame(factors[0], index = cellpair)
    interactiontensor.cc_factor = cc_df

    time_end=time.time()
    logg.info(f"Finish multiple sample SCII tensor - time cost {(((time_end - time_start) / 60) / 60)} h")

def top_pair(interactiontensor: InteractionTensor=None, pair='cc', top_n: int=15):
    """
    Identify the top CellPairs or LRpairs based on their loadings.

    Parameters
    ----------
    interactiontensor : InteractionTensor
        The InteractionTensor of the spatial transcriptomics data.
    pair : str, optional
        Type of pair to consider, either 'cc' for CellPairs or 'lr' for LRpairs, by default 'cc'.
    top_n : int, optional
        Number of top pairs to identify, by default 15.
    """
    if pair == 'cc':
        cc_df = interactiontensor.cc_factor
        top_ccpair = cc_df.loc[cc_df.sum(axis=1).sort_values(ascending=False).index[0:top_n]]
        #interactiontensor.top_ccpair=top_ccpair
        return top_ccpair
    elif pair == 'lr':
        lr_df = interactiontensor.lr_factor
        top_lrpair = lr_df.loc[lr_df.sum(axis=1).sort_values(ascending=False).index[0:top_n]]
        #interactiontensor.top_lrpair= top_lrpair
        return top_lrpair

def interaction_select_multiple(interactiontensor: InteractionTensor=None, sample: str=None, 
                                     tme_module: int=0, cellpair_module: int=0, lrpair_module: int=0, 
                                     n_lr: int=15, n_cc: int=5, 
                                ):
    """
    Plot the mean intensity heatmap for a specific TME module, CellPair module, LRPair module, and patient sample.

    Parameters
    ----------
    interactiontensor : InteractionTensor
        The InteractionTensor of the spatial transcriptomics data.
    sample : str, optional
        Patient sample identifier, by default None.
    tme_module : int, optional
        Index of the TME module to plot, by default 0.
    cellpair_module : int, optional
        Index of the CellPair module to plot, by default 0.
    lrpair_module : int, optional
        Index of the LRPair module to plot, by default 0.
    n_lr : int, optional
        Number of LRPairs to include in the plot, by default 15.
    n_cc : int, optional
        Number of CellPairs to include in the plot, by default 5.
    """

    module_lr_mt = interactiontensor.module_lr_mt
    cellpair = interactiontensor.cellpair
    lrpair = interactiontensor.lrpair
    tme_cluster = interactiontensor.tme_cluster

    individual_patient = interactiontensor.patient_map
    individual_module = interactiontensor.module_map

    factor_cc = interactiontensor.cc_factor
    factor_lr = interactiontensor.lr_factor

    top_lrpair = factor_lr.apply(lambda x: _sorted(x, n_lr))
    
    sub_cci_matrix = module_lr_mt[sample][:, :, [y for x,y,z in zip(individual_patient, individual_module, tme_cluster.tolist()) if x==sample and z == tme_module]]    
    mean_sub_mt = np.mean(sub_cci_matrix, axis=2)
    mean_df = pd.DataFrame(mean_sub_mt, index = cellpair, columns = lrpair)    
    mean_mt = mean_df[top_lrpair[lrpair_module][0:n_lr]].loc[_sorted(factor_cc[cellpair_module],n_cc).tolist()]    

    return mean_mt

def interaction_select(interactiontensor: InteractionTensor=None, 
                       tme_module: int=0, cellpair_module: int=0, lrpair_module: int=0, 
                       n_lr: int=15, n_cc: int=5, 
                      ):
    """
    Plot the mean intensity heatmap for a specific TME module, CellPair module, and LRPair module.

    Parameters
    ----------
    interactiontensor : InteractionTensor
        The InteractionTensor of the spatial transcriptomics data.
    tme_module : int, optional
        Index of the TME module to plot, by default 0.
    cellpair_module : int, optional
        Index of the CellPair module to plot, by default 0.
    lrpair_module : int, optional
        Index of the LRPair module to plot, by default 0.
    n_lr : int, optional
        Number of LRPairs to include in the plot, by default 15.
    n_cc : int, optional
        Number of CellPairs to include in the plot, by default 5.
    """

    lr_mt_list = interactiontensor.lr_mt_list_filter
    cellpair = interactiontensor.cellpair
    lrpair = interactiontensor.lrpair

    tme_cluster = interactiontensor.tme_cluster

    factor_cc = interactiontensor.cc_factor
    factor_lr = interactiontensor.lr_factor

    top_lrpair = factor_lr.apply(lambda x: _sorted(x, n_lr))
    
    sub_scii_matrix = [k for k,v in zip(lr_mt_list, pd.Series(tme_cluster).isin([tme_module])) if v]
    sub_mt = np.dstack(sub_scii_matrix)
    mean_sub_mt = np.mean(sub_mt, axis=2)
    mean_df = pd.DataFrame(mean_sub_mt, index=cellpair, columns=lrpair)
    mean_mt = mean_df[top_lrpair[lrpair_module][0:n_lr]].loc[_sorted(factor_cc[cellpair_module], n_cc).tolist()]

    return mean_mt

def merge_data(interactiontensor_list: list=None, patient_id: list=None) -> InteractionTensor:
    """
    Merge data from multiple InteractionTensor.

    Parameters
    ----------
    interactiontensor_list : list, optional
        List of InteractionTensor instances to merge, by default None.
    patient_id : list
        List of patient id

    Returns
    -------
    InteractionTensor
        A new InteractionTensor instance containing merged data.
    """
    lr_mt_list = []
    lr_mt_list_filter = []
    adata_list = []
    #zero_indices_list = []
    tme_cluster_list = []
    indice_list = []
    indice_filter_list = []
    for interactiontensor in interactiontensor_list:
        lr_mt_list.append(interactiontensor.lr_mt_list)
        lr_mt_list_filter.append(interactiontensor.lr_mt_list_filter)
        adata_list.append(interactiontensor.adata)
        #zero_indices_list.append(interactiontensor.window_zero_indices)
        tme_cluster_list.append(interactiontensor.tme_cluster)
        indice_list.append(interactiontensor.indices)
        indice_filter_list.append(interactiontensor.indices_filter)

    tmp = InteractionTensor(adata_list)

    common_columns = set(lr_mt_list_filter[0][0].columns)
    for i in range(len(lr_mt_list_filter)):
        common_columns = common_columns.intersection(set(lr_mt_list_filter[i][0].columns))
    
    common_rows = set(lr_mt_list_filter[0][0].index)
    for i in range(len(lr_mt_list_filter)):
        common_rows = common_rows.intersection(set(lr_mt_list_filter[i][0].index))
        
    tmp.lrpair = list(common_columns)
    tmp.cellpair = list(common_rows)

    for i in range(len(lr_mt_list_filter)):
        for j in range(len(lr_mt_list_filter[i])):
            lr_mt_list_filter[i][j] = lr_mt_list_filter[i][j][list(common_columns)].loc[list(common_rows)]

    cci_matrix = []
    for i in range(len(lr_mt_list_filter)):
        cci_matrix.append(np.dstack(lr_mt_list_filter[i]))

    for i in range(len(cci_matrix)):
        cci_matrix[i] = np.log1p(cci_matrix[i])

    tmp.cci_matrix_individual = cci_matrix

    indice = [item for indice in indice_filter_list for item in indice]

    tmp.indices_filter = indice

    module_lr_mt = {}
    for i in range(len(cci_matrix)):
        result = np.zeros((cci_matrix[0].shape[0], cci_matrix[0].shape[1], np.unique(tme_cluster_list[i]).max()+1))

        for cluster in range(np.unique(tme_cluster_list[i]).max()+1):
            cluster_mask = (tme_cluster_list[i] == cluster)  # Create a boolean mask for the current cluster
            cluster_data = cci_matrix[i][:, :, cluster_mask]  # Extract data for the current cluster
            result[:, :, cluster] = np.mean(cluster_data, axis=2)  # Calculate the mean along the third dimension
            module_lr_mt[patient_id[i]] = result
    tmp.module_lr_mt = module_lr_mt

    final_mt = np.nan_to_num(np.dstack(list(module_lr_mt.values())))

    tmp.cci_matrix = final_mt
 
    return tmp

def core_normalization(core, filter_threshold:bool=True, feature_range:tuple=(0, 1)):
    minv, maxv = core.min(), core.max()
    core_norm = (core-minv)/(maxv - minv)
    core_norm = core_norm*(feature_range[1] - feature_range[0]) + feature_range[0]
    if isinstance(filter_threshold, bool) and filter_threshold:
        filter_threshold=(feature_range[1]-feature_range[0])/core.shape[2]
        core_norm[core_norm < filter_threshold] = 0
    elif isinstance(filter_threshold, float):
        core_norm[core_norm < filter_threshold] = 0
    return core_norm