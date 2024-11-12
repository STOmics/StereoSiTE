import os, sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import datasets
from collections import defaultdict
from pycirclize import Circos
from pycirclize.parser import Matrix
from pycirclize.utils import calc_group_spaces, ColorCycler
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap


def _sort(s:pd.Series, num:int=100):
    sorted_s = s.sort_values(ascending=False).index[:num]
    sorted_s.index = range(num)
    return sorted_s

def scii_tensor_select(scii_tensor_list:list,
                       factor_cc:pd.DataFrame,
                       factor_lr:pd.DataFrame,
                       factor_tme:pd.DataFrame,
                       interest_TME:str, #|int
                       interest_cc_module:str, #|int
                       interest_LR_module:str, #|int
                       lr_number:int=20,
                       cc_number:int=10,
                       
                       ) -> pd.DataFrame:
    '''
    Input:
        scii_tensor_file: Name of the file that contains scii_tensor result of all windows/bins.
        factor_cc: Factor matrix of cell cell pair module, the index indicates interacitng cell pair 
                    while the column indicates cell-cell module.
        factor_lr: Factor matrix of ligand receptor pair module, the index indicates ligand receptor pair that induce cell-cell interaction
                    while the column indicates ligand-receptor module.
        factor_tme: Factor matrix of Tumor MicroEnvironment(TME) module, the index indicates window/bin
                    while the column indicates TME module.
        interest_TME: Name of the interested TME module that will be calculated.
        interest_cc_module: Name of the interested cell-cell pair module.
        interest_lr_module: Name of the interested ligand-receptor pair module.
        lr_number: The number of ligand-receptor pairs on top that will remain.
        cc_number: The number of cell-cell pair on top that will remain.
    
    Return:
        pandas.DataFrame: Matrix contains interesting interaction in the TME region of interest, index represents cell-cell pairs
                        while column represents ligand-receptor pairs.

    '''
    top_lrpair = factor_lr.apply(lambda x: _sort(x))
    tme_cluster = factor_tme.idxmax(axis=1)
    sub_scii_matrix = [k for k,v in zip(scii_tensor_list, pd.Series(tme_cluster).isin([interest_TME])) if v]
    sub_mt = np.dstack(sub_scii_matrix)
    mean_sub_mt = np.mean(sub_mt, axis=2)
    mean_df = pd.DataFrame(mean_sub_mt, index=scii_tensor_list[0].index, columns=scii_tensor_list[0].columns)
    mean_mt = mean_df[top_lrpair[interest_LR_module][0:lr_number]].loc[_sort(factor_cc[interest_cc_module], cc_number).tolist()]
    return mean_mt

def _cell_pairs_generate(cells:list, separator:str="_"):
    cell_pairs = {}
    for c1 in cells:
        for c2 in cells:
            cell_pairs[f"{c1}{separator}{c2}"] = (c1, c2)
    return cell_pairs

def scii_interaction_matrix_process(interaction_matrix:pd.DataFrame,
                                    cell_lr_separator:str="|",
                                    ) -> tuple:
    '''
    Input
        interaction_matrix: DataFrame. Matrix contains interesting interaction, index represents cell cell pairs
                            while column represents ligand receptor pairs
        cells: List of cell type names, which will be used to separate cells of the cell pair.
        separator: Separator used to combine ligand with receptor genes into LR pairs.
        cell_lr_separator: Separator used to combine cells with LR genes.
    Return
        (sectors:dict, links:list, genes:set)
        sectors: The dictionary containing all vertices and their weight. {sender-ligand: value, receiver-receptor: value, ...}
        links: List containing all links between sectors. [[sender-ligand, receiver-receptor, value], ...]
        genes: Set contains names of all genes.
    '''

    #Normalize the value into 0~100
    interaction_matrix = interaction_matrix.fillna(0).T

    scii_min, scii_max = interaction_matrix.min().min(), interaction_matrix.max().max()
    norm_mt = ((interaction_matrix - scii_min)/(scii_max - scii_min)*100).apply(round).astype(int)

    #Generate sectors and cell groups
    sectors1 = defaultdict()
    sectors2 = defaultdict()
    genes = set()
    cells = set()
    links = []
    for index, row in norm_mt.iterrows():
        sender, receiver = index
        cells.add(sender)
        cells.add(receiver)
        for LR, value in row.items():
            if value == 0:
                continue
            ligand, receptor = LR
            v1 = f"{sender}{cell_lr_separator}{ligand}"
            v2 = f"{receiver}{cell_lr_separator}{receptor}"
            if sender not in sectors1.keys():
                sectors1[sender] = defaultdict(int)
            if receiver not in sectors2.keys():
                sectors2[receiver] = defaultdict(int)
            sectors1[sender][ligand] += value
            sectors2[receiver][receptor] += value
            genes.add(ligand)
            genes.add(receptor)
            links.append([v1, v2, value])

    sectors = defaultdict(int)
    for cell in cells:
        if cell in sectors1.keys():
            for ligand, value in sectors1[cell].items():
                sectors[f"{cell}{cell_lr_separator}{ligand}"] += value
        if cell in sectors2.keys():
            for receptor, value in sectors2[cell].items():
                sectors[f"{cell}{cell_lr_separator}{receptor}"] += value
    return sectors, links, genes, list(cells)

def scii_tensor_interaction_matrix_process(interaction_matrix:pd.DataFrame,
                               cells:list,
                               separator:str="-",
                               cell_lr_separator:str="|",
                               ) -> tuple:
    '''
    Input
        interaction_matrix: DataFrame. Matrix contains interesting interaction in the TME region of interest, index represents cell-cell pairs
                            while column represents ligand-receptor pairs
        cells: list.
            List of cell type names, which will be used to separate cells of the cell pair.
        separator: str.
            Separator used to combine ligand with receptor genes into LR pairs.
        cell_lr_separator: str.
            Separator used to combine cells with LR genes.
    Return
        (sectors:dict, links:list, genes:set)
        sectors: dict.
            The dictionary containing all vertices and their weight. {sender-ligand: value, receiver-receptor: value, ...}
        links: list.
            List containing all links between sectors. [[sender-ligand, receiver-receptor, value], ...]
        genes: set.
            Set contains names of all genes.
    '''
    #Normalize the value into 0~100
    scii_min, scii_max = interaction_matrix.min().min(), interaction_matrix.max().max()
    norm_mt = ((interaction_matrix - scii_min)/(scii_max - scii_min)*100).apply(round).astype(int)

    #Generate a directory that contains names of cell types
    cell_pairs = _cell_pairs_generate(cells)

    #Generate sectors and cell groups
    sectors1 = defaultdict()
    sectors2 = defaultdict()
    genes = set()
    links = []
    for index, row in norm_mt.iterrows():
        sender, receiver = cell_pairs[index]
        for LR, value in row.items():
            if value == 0:
                continue
            ligand, receptor = LR.split(separator, 1)
            v1 = f"{sender}{cell_lr_separator}{ligand}"
            v2 = f"{receiver}{cell_lr_separator}{receptor}"
            if sender not in sectors1.keys():
                sectors1[sender] = defaultdict(int)
            if receiver not in sectors2.keys():
                sectors2[receiver] = defaultdict(int)
            sectors1[sender][ligand] += value
            sectors2[receiver][receptor] += value
            genes.add(ligand)
            genes.add(receptor)
            links.append([v1, v2, value])

    sectors = defaultdict(int)
    for cell in cells:
        if cell in sectors1.keys():
            for ligand, value in sectors1[cell].items():
                sectors[f"{cell}{cell_lr_separator}{ligand}"] += value
        if cell in sectors2.keys():
            for receptor, value in sectors2[cell].items():
                sectors[f"{cell}{cell_lr_separator}{receptor}"] += value
    return sectors, links, genes

def cells_lr_circos(interaction_matrix:pd.DataFrame,
              cells:list=None,
              cell_colors='Set3',
              gene_colors='tab20',
              link_color_palette:str='tab20',
              separater:str="-",
              cell_lr_separator:str="|",
              label_size:float=8,
              dpi:float=300,
              scii_tensor:bool=False,
              save:str=None,
              ):
    '''
    Input:
        interaction_matrix: pandas.DataFrame.
            The dataframe that contains the spatial cell interaction intensity(SCII) values of each interaction.
            The index represents cell_cell pairs and the column represents ligand_receptor pairs.
        cells: list.
            The list contains the names of all cell types. This will be used to separate the sender and receiver cells that were 
            combined to create index names for the interaction_matrix.
        cell_colors: str|dict.
            Define colors of difference cell type. User can set matplotlib's colormap (e.g. viridis, jet, tab10) or label_name -> color dict (e.g. dict(A="red", B="blue", C="green", ...)).
        gene_colors: str|dict.
            Define colors of difference gene. User can set matplotlib's colormap (e.g. viridis, jet, tab10) or label_name -> color dict (e.g. dict(A="red", B="blue", C="green", ...)).
        link_color_palette: str|list.
            Define colors of difference link. User can set matplotlib's colormap (e.g. viridis, jet, tab10) or list contain colors.
        separator: str.
            Separator used to combine ligand with receptor genes into LR pairs.
        cell_lr_separator: str.
            Separator used to combine cells with LR genes.
        scii_tensor: bool.
            If the interaction_matrix was generated by scii, set the scii_tensor=True.
        save: str.
            File name of the figure that will be saved.
    Return:    
        None
    '''

    if scii_tensor:
        sectors, links, genes = scii_tensor_interaction_matrix_process(interaction_matrix, cells=cells, separator=separater, cell_lr_separator=cell_lr_separator)
    else:
        sectors, links, genes, cells = scii_interaction_matrix_process(interaction_matrix, cell_lr_separator=cell_lr_separator)
        

    cell_groups = defaultdict(list)
    for key in sectors.keys():
        cell, gene = key.split(cell_lr_separator, 1)
        cell_groups[cell].append(key)
    group_sizes = [len(value) for key, value in cell_groups.items()]

    #Generate the color palette for cells, genes and links
    #cmap1 = plt.get_cmap('tab20b')
    #cmap2 = plt.get_cmap('tab20c')
    
    if isinstance(cell_colors, str):
        cell_color_palette = plt.get_cmap(cell_colors, len(cells)).colors
        cell_colors = dict(zip(cells, cell_color_palette[0:len(cells)]))
    elif isinstance(cell_colors, dict):
        cell_colors = cell_colors
    else:
        raise Exception("cell_colors must be string representing the matplotlib's colormap or dict")

    if isinstance(gene_colors, str):
        gene_color_palette = plt.get_cmap(gene_colors, len(genes)).colors
        gene_colors = dict(zip(sorted(list(genes)), gene_color_palette[0:len(genes)]))
    elif isinstance(gene_colors, dict):
        gene_colors = gene_colors
    else:
        raise Exception("gene_colors must be string representing the matplotlib's colormap or dict")


    lr_links = sorted(list(set([(x[0].split(cell_lr_separator, 1)[1], x[1].split(cell_lr_separator, 1)[1]) for x in links])))
    if isinstance(link_color_palette, str):
        link_color_palette = plt.get_cmap(link_color_palette, len(lr_links)).colors
    elif isinstance(link_color_palette, list):
        if len(link_color_palette) < len(lr_links):
            raise Exception("the length of link_color_palette less than then number of lr_links: {0} < {1}".format(len(link_color_palette), len(lr_links)))
    link_colors = dict(zip(lr_links, link_color_palette[0:len(lr_links)]))

    spaces = calc_group_spaces(group_sizes, space_bw_group=10, space_in_group=1)
    circos = Circos(sectors, space=spaces)

    # Plot sector track
    #ColorCycler.set_cmap("Set3")
    for sector in circos.sectors:
        track = sector.add_track(r_lim=(90, 95))
        track.axis(fc=gene_colors[sector.name.split(cell_lr_separator, 1)[1]])
        #track.text(sector.name.split("-", 1)[1], fontsize=5, r=92, orientation="vertical")

    #ColorCycler.set_cmap("tab10")
    for cell, group in cell_groups.items():
        group_deg_lim = circos.get_group_sectors_deg_lim(group)
        circos.rect(r_lim=(100, 103), deg_lim=group_deg_lim, fc=cell_colors[cell], ec="black", lw=0.5)
        group_center_deg = sum(group_deg_lim)/2
        circos.text(cell, r=106, deg=group_center_deg, adjust_rotation=True, fontsize=label_size)

    #Plot links
    for sender, receiver, value in links:
        ligand, receptor = sender.split(cell_lr_separator, 1)[1], receiver.split(cell_lr_separator, 1)[1]
        circos.link((sender, sectors[sender]-value, sectors[sender]), (receiver, sectors[receiver]-value, sectors[receiver]), color=link_colors[(ligand, receptor)], direction=1)
        sectors[sender] -= value
        sectors[receiver] -= value
        
    fig = circos.plotfig()
    #Plot legend
    rect_handles = []
    for link, color in link_colors.items():
        rect_handles.append(Patch(color=color, label=f"{link[0]}-{link[1]}"))
    rect_legend = circos.ax.legend(
        handles = rect_handles,
        bbox_to_anchor=(1.1, 1.0),
        fontsize=6,
        title="Ligand-Receptor",
    )
    circos.ax.add_artist(rect_legend)

    scatter_handles = []
    for gene, color in gene_colors.items():
        scatter_handles.append(Line2D([], [], color=color, marker="o", label=gene, ms=6, ls="None"))
    scatter_legend = circos.ax.legend(
        handles=scatter_handles,
        bbox_to_anchor=(1.5, 1.0),
        fontsize=6,
        title="Gene",
        handlelength=2,
    )
    if save != None:
        fig.savefig(save, dpi=dpi)

def cells_circos(interaction_matrix:pd.DataFrame,
                 cells:list=None,
                 cell_colors='tab20',
                 label_orientation:str="horizontal",
                 label_size:float=10,
                 dpi:float=300,
                 save:str=None,
                 scii_tensor:bool=False,
                 ):
    '''
    Input:
        interaction_matrix: pandas.DataFrame.
            The dataframe that contains the spatial cell interaction intensity(SCII) values of each interaction.
            The index represents cell_cell pairs and the column represents ligand_receptor pairs.
        cells: list.
            The list contains the names of all cell types. This will be used to separate the sender and receiver cells that were 
            combined to create index names for the interaction_matrix.
        cell_colors : str | dict[str, str], optional
            Colormap assigned to each outer track and link. User can set matplotlib's colormap (e.g. viridis, jet, tab10) or label_name -> color dict (e.g. dict(A="red", B="blue", C="green", ...))
        save: str.
            File name of the figure that will be saved.
        scii_tensor: bool.
            If the interaction_matrix was generated by scii, set the scii_tensor=True.
    Return:
        None                       
    '''
    #Normalize the value into 0~100
    scii_min, scii_max = interaction_matrix.min().min(), interaction_matrix.max().max()
    norm_interaction_matrix = ((interaction_matrix - scii_min)/(scii_max - scii_min)*100).apply(round).astype(int)

    #Generate the matrix that will be used to draw circos
    if scii_tensor:
        if not isinstance(cells, list):
            cells = [cell for cell in cells]
        if cells==None:
            raise Exception("When scii_tensor is True, the cells parameter must be given")
        cci_df = norm_interaction_matrix.sum(axis=1).to_frame()
        cell_pairs = _cell_pairs_generate(cells)
        cci_df[['sender', 'receiver']] = [cell_pairs[x] for x in cci_df.index]
        
    else:
        cci_df = norm_interaction_matrix.T.sum(axis=1).to_frame()
        cci_df[['sender', 'receiver']] = [list(x) for x in cci_df.index]    
    
    cells = list(set(cci_df['sender'].unique()) | set(cci_df['receiver'].unique()))
    cells_dict = dict(zip(cells, range(len(cells))))
    cci_df['sender_index'] = cci_df['sender'].map(cells_dict)
    cci_df['receiver_index'] = cci_df['receiver'].map(cells_dict)

    fromto_table_df = cci_df[['sender', 'receiver', 0]].rename(columns = {'sender': 'from', 'receiver': 'to', 0: 'value'}).reset_index(drop=True)
    matrix = Matrix.parse_fromto_table(fromto_table_df)

    #Draw circos
    circos = Circos.initialize_from_matrix(
        matrix,
        space=3,
        cmap=cell_colors,
        #ticks_interval=5,
        label_kws=dict(size=label_size, r=110, orientation=label_orientation),
        link_kws=dict(direction=1, ec='black', lw=0.5),
    )
    fig = circos.plotfig()
    if save != None:
        fig.savefig(save, dpi=dpi)
    
def lr_circos(interaction_matrix:pd.DataFrame,
              cells:list,
              cmap='Set3',
              separator:str="-",
              label_orientation:str="vertical",
              label_size:float=6,
              dpi:float=300,
              save:str=None,
              scii_tensor=False,
              ):
    '''
    Input:
        interaction_matrix: The dataframe that contains the spatial cell interaction intensity(SCII) values of each interaction.
                            The index represents cell_cell pairs and the column represents ligand_receptor pairs.
        cmap : str | dict[str, str], optional
            Colormap assigned to each outer track and link. User can set matplotlib's colormap (e.g. viridis, jet, tab10) or label_name -> color dict (e.g. dict(A="red", B="blue", C="green", ...))
        separator: str.
            Separator used to combine ligand with receptor genes into LR pairs.
        save: File name of the figure that will be saved.
    '''
    if scii_tensor:
        _, links, _ = scii_tensor_interaction_matrix_process(interaction_matrix, cells, separator=separator)
    else:
        _, links, _ = scii_interaction_matrix_process(interaction_matrix)
    links_df = pd.DataFrame(links, columns = ['from', 'to', 'value']).sort_values(by = ['from', 'to'])
    matrix = Matrix.parse_fromto_table(links_df)
    circos = Circos.initialize_from_matrix(
        matrix,
        space=2,
        cmap=cmap,
        #ticks_interval=5,
        label_kws=dict(size=label_size, r=110, orientation=label_orientation),
        link_kws=dict(direction=1, ec='black', lw=0.5, alpha=0.5),
    )
    fig = circos.plotfig()
    if save != None:
        fig.savefig(save, dpi=dpi)



