import os, sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import cm

def core_heatmap(core,
                save:str=None,
                figsize:tuple=(10, 6),
                dpi=300,
                cmap:str='rocket',
                ):
    fig = plt.figure(111, figsize=figsize, facecolor='white')
    ax1 = fig.add_subplot(1, 2, 1, projection='3d')
    ax11 = fig.add_subplot(1, 2, 2) #color bar

    # Set the azimuth
    ax1.azim = -40
    ax1.dist = 9
    ax1.elev = 25

    # Draw grid
    x = np.arange(core.shape[0])
    y = np.arange(core.shape[1])
    X, Y = np.meshgrid(x, y)

    # Draw the 2D heatmap
    for module in range(core.shape[2]):
        paper_white = cm.ScalarMappable(cmap = cmap).to_rgba(core[:,:,module])
        surf = ax1.plot_surface(X, np.zeros(shape=X.shape)+module, Y, rstride=1, cstride=1, facecolors=paper_white, linewidth=0, antialiased=True, alpha=1)

    # Draw main coordinates and set the ticks
    ax1.tick_params(axis='x', colors='k', labelsize=6)
    ax1.tick_params(axis='y', colors='k', labelsize=6)
    ax1.tick_params(axis='z', colors='k', labelsize=6)
    ax1.set_xlabel("LR Modules")
    ax1.set_ylabel("TME Modules")
    ax1.set_zlabel("CC Modules")
    ax1.set_xticks(range(core.shape[0]))
    ax1.set_yticks(range(core.shape[2]))
    ax1.set_zticks(range(core.shape[1]))

    # Set the color bar
    ax11.set_visible(False)
    axin11 = inset_axes(ax11, width="2%", height="75%", loc='center left', borderpad=0)
    axin11.tick_params(axis='y', labelsize=12)
    norm_cot = mpl.colors.Normalize(vmin=core.min(), vmax=core.max())
    fig.colorbar(cm.ScalarMappable(norm=norm_cot,cmap=cmap), shrink=1.0, aspect=5, ax=ax11, cax=axin11)
    if save != None:
        plt.savefig(save, dpi=dpi)
    else:
        plt.show()
        plt.close()

def tme_core_heatmap(core,
                     tme_number:int,
                     vmin:int=0,
                     vmax:int=None,
                     save:str=None,
                     figsize:tuple=(5, 4),
                     dpi:float=300,
                     cmap:str='rocket',
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
                    cmap=cmap,
                     )
    cbar = h.figure.colorbar(h.collections[0])
    ticks = [tme_df.min().min(), vmax] # tme7_df.max().max()]
    labels = ['Low', 'High']
    cbar.set_ticks(ticks)
    cbar.set_ticklabels(labels)
    plt.title(f"TME Module {tme_number}")
    if save != None:
        plt.savefig(save, dpi=dpi)
    else:
        plt.show()
        plt.close()

def interaction_heatmap(interaction_matrix:pd.DataFrame,
                        linewidths:float=0.005,
                        linecolor:str='black',
                        vmax:int=None,
                        save:str=None,
                        figsize:tuple=(6, 1.5),
                        dpi:float=300,
                        cmap:str='rocket',
                        ):
    if vmax==None:
        vmax=interaction_matrix.max().max()
    plt.figure(figsize=figsize)
    h = sns.heatmap(interaction_matrix,
                    linewidths=linewidths,
                    linecolor=linecolor,
                    vmax=vmax,
                    cbar=False,
                    cmap=cmap,
                    )
    cbar = h.figure.colorbar(h.collections[0])
    ticks = [0, vmax]
    labels = ['Low', 'High']
    cbar.set_ticks(ticks)
    cbar.set_ticklabels(labels)
    if save != None:
        plt.savefig(save, dpi=dpi)
    else:
        plt.show()
        plt.close()

def reconstruction_error_line(re_mat:np.ndarray,
                              figsize:tuple=(4, 4),
                              save:str=None,
                              palette:str='tab20',
                              dpi=300,
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
        plt.savefig(save, bbox_inches='tight', dpi=dpi)
    else:
        plt.show()
        plt.close()