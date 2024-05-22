import os, sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

def tme_core_heatmap(core,
                     tme_number:int,
                     vmin:int=0,
                     vmax:int=None,
                     save:str=None,
                     figsize:tuple=(5, 4),
                     ):
    tme_df = pd.DataFrame(core[:, :, tme_number])
    tme_df.columns = tme_df.columns.map(lambda x: f"LR_Module {x}")
    tme_df.index = tme_df.index.map(lambda x: f"CC_Module {x}")
    if vmax==None:
        vmax=tme_df.max().max()
    plt.figure(figsize=figsize)
    h = sns.heatmap(tme_df, #cmap='Purples',
                    linewidths=0.005, linecolor='black',
                    annot=False,
                    cbar=False,
                    vmin=vmin,
                    vmax=vmax,
                     )
    cbar = h.figure.colorbar(h.collections[0])
    ticks = [tme_df.min().min(), vmax] # tme7_df.max().max()]
    labels = ['Low', 'High']
    cbar.set_ticks(ticks)
    cbar.set_ticklabels(labels)
    plt.title(f"TME Module {tme_number}")
    if save != None:
        plt.savefig(save)
    else:
        plt.show()
        plt.close()

def interaction_heatmap(interaction_matrix:pd.DataFrame,
                        linewidths:float=0.005,
                        linecolor:str='black',
                        vmax:int=None,
                        save:str=None,
                        figsize:tuple=(6, 1.5),
                        ):
    if vmax==None:
        vmax=interaction_matrix.max().max()
    plt.figure(figsize=figsize)
    h = sns.heatmap(interaction_matrix,
                    linewidths=linewidths,
                    linecolor=linecolor,
                    vmax=vmax,
                    cbar=False
                    )
    cbar = h.figure.colorbar(h.collections[0])
    ticks = [0, vmax]
    labels = ['Low', 'High']
    cbar.set_ticks(ticks)
    cbar.set_ticklabels(labels)
    if save != None:
        plt.savefig(save)
    else:
        plt.show()
        plt.close()

def reconstruction_error_line(re_mat:np.ndarray,
                              figsize:tuple=(4, 4),
                              save:str=None,
                              palette:str='tab20',
                              ):

    num_TME_modules = re_mat.shape[0]
    colors = plt.get_cmap(palette, num_TME_modules-2).colors
 
    plt.figure(figsize=figsize)
    for j in range(3,num_TME_modules):
        plt.plot(np.arange(3, len(re_mat[j])),re_mat[j][3:],label = 'rank: TME={}'.format(j), color=colors[j-2])
    plt.xlabel('ranks: CC, LR')
    plt.ylabel('reconstruction error')
    plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left")
    if save!=None:
        plt.savefig(save, bbox_inches='tight')
    else:
        plt.show()
        plt.close()