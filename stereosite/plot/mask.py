#! /usr/bin/env python3
# _*_ coding: utf-8 _*_

#@Author:  LiuXing liuxing2@genomics.cn
#@Date:  2023-07-11 16:04:39
#@Last Modified by: LiuXing
#@Last Modified time: 2023-07-11 16:04:39

import numpy as np
import tifffile
import cv2
import os, sys
import anndata

def _color_change(value):
    digit = list(map(str, range(10))) + list("abcdef")
    a1 = digit.index(value[1]) * 16 + digit.index(value[2])
    a2 = digit.index(value[3]) * 16 + digit.index(value[4])
    a3 = digit.index(value[5]) * 16 + digit.index(value[6])
    return (a3, a2, a1)

def _get_color_list(my_color):
    if (my_color == None):
        my_color = ['#0343df', '#f97306', '#15b01a', '#e50000', '#7e1e9c', '#653700', '#ff81c0', '#929591', 
            '#6e750e', '#00ffff', '#ff796c', '#06c2ac', '#75bbfd', '#01ff07', '#cb416b', '#bf77f6', '#ceb301', 
            '#137e6d', '#516572', '#dbb40c', '#d0fefe', '#9ffeb0', '#fdaa48', '#ffcfdc', '#ffffc2', '#ac9362', 
            '#7a9703', '#96ae8d', '#b66a50', '#411900']
    rgb_mycolor = []
    for color in my_color:
        rgb_mycolor.append(list(_color_change(color)))
    return rgb_mycolor

def mask_coloring(adata: anndata, mask_file: str, 
                  anno: str='cell2loc_anno', 
                  save: str=None, save_legend: str=None):
    """
    Paint each cell in the mask with color corresponding to the palette of cell type in anndata

    Parameters
    -----------
    adata
        anndata file with cell type annotation
    mask_file
        cell mask file path
    anno
        annotation key in adata. default=cell2loc_anno
    save
        if save was specified with a file path, the colored mask will be writen to it. The file suffix should be .png|.jpg
    save_legend
        if save_legend was specified with a file path, the corresponding legend of colored mask will be writen to it. The file suffix should be .png|.jpg

    Returns
    ----------
    Tuple <colored_cell_mask, legend>
        colored cell mask: ndarray
        legend: ndarray
    """
    maskImg = tifffile.imread(mask_file)
    if (maskImg.max() == 1):
        _, labels = cv2.connectedComponents(maskImg)
    else:
        labels = maskImg
    dst = np.nonzero(labels)
    paletteKey = f"{anno}_colors"
    my_color = list(adata.uns[paletteKey]) if paletteKey in adata.uns.keys() else None
    rgb_color = _get_color_list(my_color)
    clusterDict = dict(zip(adata.obs.index.astype(int), adata.obs[anno].cat.codes))
    img_colors = [rgb_color[int(clusterDict[x])] if x in clusterDict.keys() else [255, 255, 255] for x in labels[dst]]
    new_img = np.zeros((maskImg.shape[0], maskImg.shape[1], 3), dtype = np.int8)
    new_img[dst] = img_colors

    #get legend
    cells = adata.obs[anno].cat.categories
    width = 10000
    high = 1000
    legend = np.zeros([high*len(cells), width,3])
    legend.fill(255)
    font = cv2.FONT_HERSHEY_SIMPLEX
    fontScale = 20  
    # Line thickness of 40 px
    thickness = 50
    pointSize = 100
    for i in range(len(cells)):
        color = rgb_color[i]
        cell = cells[i]
        point = (int(width/8), int((i+0.4)*high))
        cv2.circle(legend, point, pointSize, color, pointSize*2)
        org = (int(width/5), int((i+0.6)*high))
        legend = cv2.putText(legend, cell, org, font, fontScale, (0, 0, 0), thickness, cv2.LINE_AA)
    
    if save != None:
        out_dir = os.path.dirname(save)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir, exist_ok=True)
        cv2.imwrite(save, new_img)
    if save_legend != None:
        out_dir = os.path.dirname(save_legend)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        cv2.imwrite(save_legend, legend)
    return new_img, legend