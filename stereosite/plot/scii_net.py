import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import igraph as ig
from igraph import Graph
from .scii_circos import interaction_matrix_process

def lr_link_graph_generate(interaction_matrix:pd.DataFrame,
                  cells:list,
                  separator:str="-",
                  cell_lr_separator="|",
                  reducer:int=3,
                  ) -> Graph:
    '''
    Input:
        interaction_matrix: DataFrame. Matrix contains interesting interaction in the TME region of interest, index represents cell-cell pairs
                            while column represents ligand-receptor pairs
        cells: List of cell type names, which will be used to separate cells of the cell pair.
        separator: Separator used to combine ligand with receptor genes into LR pairs.
        cell_lr_separator: Separator used to combine cells with LR genes.
        reducer: The size of a vertex = weight of the vertex/reducer
    Return:
        g: Graph containing all vertices(ligand or receptor) and their links information.
        
    '''
    #generate color palette for cells and links
    cell_palette = 'Set3'
    cell_colors = dict(zip(cells, [tuple(x[0:3]) for x in plt.get_cmap(cell_palette, len(cells)).colors]))

    #transfer matrix into igraph
    sectors, links, genes = interaction_matrix_process(interaction_matrix, cells, separator=separator, cell_lr_separator=cell_lr_separator)
    vertices = list(sectors.keys())
    vertices_color = [cell_colors[x.split(cell_lr_separator)[0]] for x in vertices]
    vertices_dict = dict(zip(vertices, range(len(vertices))))
    edges_index = [(vertices_dict[x[0]], vertices_dict[x[1]]) for x in links]
    lr_links = set([(x[0].split(cell_lr_separator)[1], x[1].split(cell_lr_separator)[1]) for x in links])
    lr_palette = "tab20"
    link_color_palette = [tuple(x[0:3]) for x in plt.get_cmap(lr_palette, len(lr_links)).colors]
    link_colors = dict(zip(lr_links, link_color_palette))
    edges_color = [link_colors[(x[0].split(cell_lr_separator)[1], x[1].split(cell_lr_separator)[1])] for x in links]

    g = Graph(n=len(vertices), edges=edges_index, directed=True)
    g.vs['name'] = vertices
    g.vs['label'] = [x.split(cell_lr_separator, 1)[1] for x in vertices]
    g.vs['color'] = vertices_color
    g.vs['weight'] = [x/reducer if x>reducer*5 else 5 for x in sectors.values()]
    g.es['weight'] = [x[2] for x in links]
    g.es['color'] = edges_color
    return g

def grap_plot(interaction_matrix:pd.DataFrame,
              cells:list,
              separator:str='-',
              cell_lr_separator:str="|",
              layout_type:str='kk',
              save:str=None,
              **kwargs,
              ):
    '''
    Input:
        interaction_matrix: DataFrame. Matrix contains interesting interaction in the TME region of interest, index represents cell-cell pairs
                            while column represents ligand-receptor pairs
        cells: List of cell type names, which will be used to separate cells of the cell pair.
        separator: Separator used to combine ligand with receptor genes into LR pairs.
        cell_lr_separator: Separator used to combine cells with LR genes.
        layout_type: Layout style used to draw the graph. Same to the layout of igraph
        save: File name of the figure that will be saved.
        **kwargs: Dictionary containing parameters of igraph plot setting.
    Return:
        None
    '''
    
    g = lr_link_graph_generate(interaction_matrix, cells, separator=separator, cell_lr_separator=cell_lr_separator)
    cell_colors = dict(zip([x.split(cell_lr_separator)[0] for x in g.vs['name']], g.vs['color']))
    link_colors = dict(zip([f"{x.source_vertex['label']}{separator}{x.target_vertex['label']}" for x in g.es], g.es['color']))

    fig, ax = plt.subplots()
    layout = g.layout(layout_type)
    layout.rotate(0)
    #draw graph
    ig.plot(g,
            layout = layout,
            vertex_size=g.vs['weight'],
            vertex_label_angle=90,
            vertex_label_dist = 0,
            vertex_label_size = 8,
            edge_width=[0.5, 4],
            edge_curved=0.2,
            edge_arrow_size=10,
            edge_arrow_width=5,
            target=ax,
            )
    #generate legend
    cell_legend_handles = []
    for cell, color in cell_colors.items():
        handle = ax.scatter(
            [], [],
            s=100,
            facecolor=color,
            label=cell,
        )
        cell_legend_handles.append(handle)
    l1 = ax.legend(
        handles=cell_legend_handles,
        title='Cell Type',
        bbox_to_anchor=(1.0, 1.0),
        bbox_transform=ax.transAxes,
    )
    lr_legend_handles = []
    for lr, color in link_colors.items():
        handle = ax.scatter(
            [], [],
            s=100,
            facecolor=color,
            label=lr,
            marker = 's',
        )
        lr_legend_handles.append(handle)
    ax.legend(
        handles=lr_legend_handles,
        title='Ligand-Receptor',
        bbox_to_anchor=(1.0, 0.6),
        bbox_transform=ax.transAxes,
    )
    fig.gca().add_artist(l1)
    fig.set_size_inches(15, 15)
    if save != None:
        fig.savefig(save)
