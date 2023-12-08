import os, sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import datasets
from collections import defaultdict

class PySankeyException(Exception):
    pass


class NullsInFrame(PySankeyException):
    pass


class LabelMismatch(PySankeyException):
    pass

def _check_data_matches_labels(labels, data):
    if isinstance(data, list):
        data = set(data)
    if isinstance(labels, pd.Series):
        data = set(data.unique().tolist())
    if isinstance(labels, list):
        labels = set(labels)
    if labels != data:
        raise LabelMismatch('Labels and data do not match. {0}'.format(','.join(labels)))

def core_process(core:np.array, model_names:list = ['CC_Module', 'LR_Module', 'TME']):
    # preprocess core matrix with three dimensionality: ligand-receptor model, cell-cell model, TME model
    left = []
    mid_l = []
    mid_r = []
    right = []
    left_weight = []
    mid_l_weight = []
    mid_r_weight = []
    right_weight = []
    lr_tme_link = core.sum(axis=0)
    cc_tme_link = core.sum(axis=1)
    lrs, ccs, tmes = core.shape
    for tme in range(tmes):
        for lr in range(lrs):
            left.append(f"{model_names[0]} {lr}")
            mid_l.append(f"{model_names[2]} {tme}")
            left_weight.append(cc_tme_link[lr, tme])
            #mid_l_weight.append(cc_tme_link[lr, tme])
        for cc in range(ccs):
            right.append(f"{model_names[1]} {cc}")
            mid_r.append(f"{model_names[2]} {tme}")
            right_weight.append(lr_tme_link[cc, tme])
            #mid_r_weight.append(lr_tme_link[cc, tme])
    mid_l_weight = tmes*ccs*[cc_tme_link.sum()/(tmes*ccs)]
    mid_r_weight = tmes*lrs*[lr_tme_link.sum()/(tmes*lrs)]
    
    dataFrame = pd.DataFrame({'left': left, 'right': right, 'mid_l': mid_l, 'mid_r': mid_r, 
                          'left_weight': left_weight, 'right_weight': right_weight, 
                              'mid_l_weight': mid_l_weight, 'mid_r_weight': mid_r_weight}, index=range(len(left)))
    return dataFrame

def factor_process(factors:np.array, names:list, model_type:str='CC_Module'):
    left = []
    right = []
    leftWeight = []
    rightWeight = []
    for index, element in np.ndenumerate(factors):
        left.append(f'{model_type} {index[1]}')
        right.append(names[index[0]])
        leftWeight.append(element)
        rightWeight.append(element)
    dataFrame = pd.DataFrame({'left': left, 'right': right, 'leftWeight': leftWeight, 'rightWeight': rightWeight}, index=range(len(left)))
    return dataFrame

def sankey_3d(data:pd.DataFrame, color_dict=None, left_labels = None, right_labels = None, mid_labels = None,
              aspect=3, fontsize=5, save=None, close_plot=False, patch_alpha:float=0.99, link_alpha:float=0.4,
              palette:str='tab20', interval:float=0.005):
    '''
    Make Sankey Diagram showing flow:  ligand-receptor model <--- TME ---> cell-cell model
    
    Inputs:
        data: pandas.dataFrame contains columns left, right, mid_l, mid_r, left_weight, right_weight, mid_l_weight, mid_r_weight
        color_dict: Dictionary of colors to use for each label {'label': 'color'}
        left_labels: Order of the left labels in the diagram
        right_labels: Order of the right labels in the diagram
        mid_labels: Order of the middle labels in the diagram
        aspect: Vertical extent of the diagram in units of horizontal extent
        fontsize: Fontsize of patch label text
        save: If the figure file name was given, the sankey figure will be stored in it.
        palette: Palette style
        interval: Distance between two adjacent patchs = interval * vertical length of all patchs.
        
    Output:
        None
    '''
    plt.figure(dpi=140)
    plt.rc('text', usetex=False)
    plt.rc('font', family='serif')

    if len(data[(data.left.isnull()) | (data.right.isnull()) | (data.mid_l.isnull()) | (data.mid_r.isnull())]):
        raise NullsInFrame('Sankey graph dose not support null values.')

    #Identify all labels that appear 'left' or 'right'
    allLabels = pd.Series(np.r_[data.mid_r.unique(), data.left.unique(), data.right.unique(), data.mid_l.unique()]).unique()

    #Identify labels
    if left_labels == None:
        left_labels = data.left.unique()[::-1]
    else:
        _check_data_matches_labels(left_labels, data['left'])
    if mid_labels == None:
        mid_labels = data.mid_r.unique()[::-1]
    else:
        _check_data_matches_labels(mid_labels, data['mid_l'])
    if right_labels == None:
        right_labels = data.right.unique()[::-1]
    else:
        _check_data_matches_labels(right_labels, data['right'])

    if color_dict is None:
        color_dict = {}
        colorPalette = sns.color_palette(palette, len(allLabels))
        for i, label in enumerate(allLabels):
            color_dict[label] = colorPalette[i]
        #for label in left_labels:
        #    color_dict[label] = colorPalette[i+1]
        #for label in right_labels:
        #    color_dict[label] = colorPalette[i+2]
    else:
        missing = [label for label in allLabels if label not in color_dict.keys()]
        if missing:
            msg = "The color_dict parameter is missing values for the following labels: {}".format(', '.join(missing))
            raise ValueError(msg)

    #Determine widths of individual strips
    
    ns_l = defaultdict()
    ns_m_l = defaultdict()
    ns_m_r = defaultdict()
    ns_r = defaultdict()
    for midLabel in mid_labels:
        leftDict = {}
        midLDict = {}
        midRDict = {}
        rightDict = {}
        for leftLabel in left_labels:
            midLDict[leftLabel] = data[(data.mid_l == midLabel) & (data.left == leftLabel)].mid_l_weight.sum()
            leftDict[leftLabel] = data[(data.mid_l == midLabel) & (data.left == leftLabel)].left_weight.sum()
        ns_m_l[midLabel] = midLDict
        ns_l[midLabel] = leftDict
        for rightLabel in right_labels:
            midRDict[rightLabel] = data[(data.mid_r == midLabel) & (data.right == rightLabel)].mid_r_weight.sum()
            rightDict[rightLabel] = data[(data.mid_r == midLabel) & (data.right == rightLabel)].right_weight.sum()
        ns_m_r[midLabel] = midRDict
        ns_r[midLabel] = rightDict

    midLWidths = defaultdict()
    for i, midLabel in enumerate(mid_labels):
        myD = {}
        myD['mid'] = data[data.mid_l == midLabel].mid_l_weight.sum()
        if i == 0:
            myD['bottom'] = 0
            myD['top'] = myD['mid']
        else:
            myD['bottom'] = midLWidths[mid_labels[i - 1]]['top'] + interval*data.mid_l_weight.sum()
            myD['top'] = myD['bottom'] + myD['mid']
            topEdge = myD['top']
        midLWidths[midLabel] = myD
    midRWidths = defaultdict()
    for i, midLabel in enumerate(mid_labels):
        myD = {}
        myD['mid'] = data[data.mid_r == midLabel].mid_r_weight.sum()
        if i == 0:
            myD['bottom'] = 0
            myD['top'] = myD['mid']
        else:
            myD['bottom'] = midRWidths[mid_labels[i - 1]]['top'] + interval*data.mid_r_weight.sum()
            myD['top'] = myD['bottom'] + myD['mid']
            topEdge = myD['top']
        midRWidths[midLabel] = myD

    # Determine positions of left label patches and total widths
    leftWidths = defaultdict()
    for i, leftLabel in enumerate(left_labels):
        myD = {}
        myD['left'] = data[data.left == leftLabel].left_weight.sum()
        if i == 0:
            myD['bottom'] = 0
            myD['top'] = myD['left']
        else:
            myD['bottom'] = leftWidths[left_labels[i - 1]]['top'] + interval*data.left_weight.sum()
            myD['top'] = myD['bottom'] + myD['left']
            topEdge = myD['top']
        leftWidths[leftLabel] = myD

    # Determine positions of right label patches and total widths
    rightWidths = defaultdict()
    for i, rightLabel in enumerate(right_labels):
        myD = {}
        myD['right'] = data[data.right == rightLabel].right_weight.sum()
        if i == 0:
            myD['bottom'] = 0
            myD['top'] = myD['right']
        else:
            myD['bottom'] = rightWidths[right_labels[i-1]]['top'] + interval*data.right_weight.sum()
            myD['top'] = myD['bottom'] + myD['right']
            topEdge = myD['top']
        rightWidths[rightLabel] = myD

    #Total vertical extent of diagram
  
    xMax = topEdge/aspect
    width = 0.35*xMax
    pad = 0.9995

    #Draw vertical bars on left and right of each label's section & print label
    for midLabel in mid_labels:
        plt.fill_between(
            [-width/2, width/2],
            2 * [midLWidths[midLabel]['bottom'] *pad],
            2 * [(midLWidths[midLabel]['bottom'] + midLWidths[midLabel]['mid']) * pad],
            color = color_dict[midLabel],
            alpha=patch_alpha
        )
        plt.text(
            0,
            midLWidths[midLabel]['bottom'] + 0.5*midLWidths[midLabel]['mid'],
            midLabel,
            {'ha': 'center', 'va': 'center'},
            fontsize = fontsize
        )

    for leftLabel in left_labels:
        plt.fill_between(
            [-xMax-width, -xMax],
            2 * [leftWidths[leftLabel]['bottom'] * pad],
            2 * [(leftWidths[leftLabel]['bottom'] + leftWidths[leftLabel]['left']) * pad],
            color = color_dict[leftLabel],
            alpha=patch_alpha
        )
        plt.text(
            -xMax - width/2,
            leftWidths[leftLabel]['bottom'] + 0.5*leftWidths[leftLabel]['left'],
            leftLabel,
            {'ha': 'center', 'va': 'center'},
            fontsize=fontsize
        )

    for rightLabel in right_labels:
        plt.fill_between(
            [xMax, xMax + width],
            2 * [rightWidths[rightLabel]['bottom'] * pad],
            2 * [(rightWidths[rightLabel]['bottom'] + rightWidths[rightLabel]['right']) * pad],
            color = color_dict[rightLabel],
            alpha=patch_alpha
        )
        plt.text(
            xMax + width/2,
            rightWidths[rightLabel]['bottom'] + 0.5*rightWidths[rightLabel]['right'],
            rightLabel,
            {'ha': 'center', 'va': 'center'},
            fontsize = fontsize
        )

    # Plot strips
    for midLabel in mid_labels:
        for leftLabel in left_labels:
            labelColor = midLabel
            if len(data[(data.mid_l == midLabel)& (data.left == leftLabel)]) > 0:
                # Create array of y values for each strip, half at middle value, half at left
                ys_d = np.array(50 * [leftWidths[leftLabel]['bottom']] + 50*[midLWidths[midLabel]['bottom']])
                ys_d = np.convolve(ys_d, 0.05*np.ones(20), mode='valid')
                ys_d = np.convolve(ys_d, 0.05*np.ones(20), mode='valid')
                ys_u = np.array(50*[leftWidths[leftLabel]['bottom'] + ns_l[midLabel][leftLabel]] + 50*[midLWidths[midLabel]['bottom'] + ns_m_l[midLabel][leftLabel]])
                ys_u = np.convolve(ys_u, 0.05 * np.ones(20), mode='valid')
                ys_u = np.convolve(ys_u, 0.05 * np.ones(20), mode='valid')

                # Update bottom edges at each label so next strip starts at the right place
                midLWidths[midLabel]['bottom'] += ns_m_l[midLabel][leftLabel]
                leftWidths[leftLabel]['bottom'] += ns_l[midLabel][leftLabel]
                plt.fill_between(
                    np.linspace(-xMax, -width/2, len(ys_d)), ys_d, ys_u, alpha=link_alpha,
                    color = color_dict[labelColor]
                )

        for rightLabel in right_labels:
            labelColor = midLabel
            if len(data[(data.mid_r == midLabel) & (data.right == rightLabel)]) > 0 :
                # Create array of y values for each strip, half at let value,
                # half at right
                ys_d = np.array(50 * [midRWidths[midLabel]['bottom']] + 50 * [rightWidths[rightLabel]['bottom']])
                ys_d = np.convolve(ys_d, 0.05*np.ones(20), mode='valid')
                ys_d = np.convolve(ys_d, 0.05*np.ones(20), mode='valid')
                ys_u = np.array(50 * [midRWidths[midLabel]['bottom'] + ns_m_r[midLabel][rightLabel]] + 50 * [rightWidths[rightLabel]['bottom'] + ns_r[midLabel][rightLabel]])
                ys_u = np.convolve(ys_u, 0.05 * np.ones(20), mode='valid')
                ys_u = np.convolve(ys_u, 0.05 * np.ones(20), mode='valid')

                # Update bottom edges at each label so next strip starts at the right place
                midRWidths[midLabel]['bottom'] += ns_m_r[midLabel][rightLabel]
                rightWidths[rightLabel]['bottom'] += ns_r[midLabel][rightLabel]
                plt.fill_between(
                    np.linspace(width/2, xMax, len(ys_d)), ys_d, ys_u, alpha=link_alpha,
                    color = color_dict[labelColor]
                )
    plt.gca().axis('off')
    plt.gcf().set_size_inches(6, 6)
    if save != None:
        plt.savefig(save, bbox_inches='tight', dpi=150)
    if close_plot:
        plt.close()

def sankey_2d(data:pd.DataFrame, color_dict=None, left_labels = None, right_labels = None,
              aspect=4, fontsize=4, save=None, close_plot=False, patch_alpha:float=0.99, link_alpha:float=0.65,
              interval:float=0.0
              ):
    '''
    Make Sankey Diagram
    
    Inputs:
        data: pandas.DataFrame contains columns left, right, mid_l, mid_r, left_weight, right_weight, mid_l_weight, mid_r_weight
        color_dict: Dictionary of colors to use for each label {'label': 'color'}
        left_labels: order of the left labels in the diagram
        right_labels: order of the right labels in the diagram
        aspect: vertical extent of the diagram in units of horizontal extent
        
    Output:
        None
    '''  
    plt.figure(dpi=140)
    plt.rc('text', usetex=False)
    plt.rc('font', family='serif')

    if len(data[(data.left.isnull()) | (data.right.isnull())]):
        raise NullsInFrame('Sankey graph dose not support null values.')

    #Identify all labels that appear 'left' or 'right'
    allLabels = pd.Series(np.r_[data.left.unique(), data.right.unique()]).unique()

    #Identify left labels
    if left_labels == None:
        left_labels = data.left.unique()[::-1]
    else:
        _check_data_matches_labels(left_labels, data['left'])
    if right_labels == None:
        right_labels = data.right.unique()[::-1]
    else:
        _check_data_matches_labels(right_labels, data['right'])

    color_dict = {}
    palette = 'tab20'
    colorPalette = sns.color_palette(palette, len(allLabels))
    for i, label in enumerate(allLabels):
        color_dict[label] = colorPalette[i]

    #Determine widths of individual strips
    from collections import defaultdict
    ns_l = defaultdict()
    ns_r = defaultdict()
    for leftLabel in left_labels:
        leftDict = {}
        rightDict = {}
        for rightLabel in right_labels:
            leftDict[rightLabel] = data[(data.left == leftLabel) & (data.right == rightLabel)].leftWeight.sum()
            rightDict[rightLabel] = data[(data.left == leftLabel) & (data.right == rightLabel)].rightWeight.sum()
        ns_l[leftLabel] = leftDict
        ns_r[leftLabel] = rightDict

    # Determine positions of left label patches and total widths
    leftWidths = defaultdict()
    for i, leftLabel in enumerate(left_labels):
        myD = {}
        myD['left'] = data[data.left == leftLabel].leftWeight.sum()
        if i == 0:
            myD['bottom'] = 0
            myD['top'] = myD['left']
        else:
            myD['bottom'] = leftWidths[left_labels[i - 1]]['top'] + interval*data.leftWeight.sum()
            myD['top'] = myD['bottom'] + myD['left']
            topEdge = myD['top']
        leftWidths[leftLabel] = myD

    # Determine positions of right label patches and total widths
    rightWidths = defaultdict()
    for i, rightLabel in enumerate(right_labels):
        myD = {}
        myD['right'] = data[data.right == rightLabel].rightWeight.sum()
        if i == 0:
            myD['bottom'] = 0
            myD['top'] = myD['right']
        else:
            myD['bottom'] = rightWidths[right_labels[i-1]]['top'] + interval*data.rightWeight.sum()
            myD['top'] = myD['bottom'] + myD['right']
            topEdge = myD['top']
        rightWidths[rightLabel] = myD

    #Total vertical extent of diagram
    l_width = -0.3
    r_width = 1.3
    xMax = topEdge/aspect

    #Draw vertical bars on left and right of each label's section & print label
    for leftLabel in left_labels:
        plt.fill_between(
            [l_width*xMax, 0],
            2 * [leftWidths[leftLabel]['bottom']],
            2 * [leftWidths[leftLabel]['bottom'] + leftWidths[leftLabel]['left']],
            color = color_dict[leftLabel],
            alpha=patch_alpha
        )
        plt.text(
            l_width/2 * xMax,
            leftWidths[leftLabel]['bottom'] + 0.5*leftWidths[leftLabel]['left'],
            leftLabel,
            {'ha': 'center', 'va': 'center'},
            fontsize=fontsize
        )
    for rightLabel in right_labels:
        plt.fill_between(
            [xMax, r_width*xMax], 2 * [rightWidths[rightLabel]['bottom']],
            2 * [rightWidths[rightLabel]['bottom'] + rightWidths[rightLabel]['right']],
            color = color_dict[rightLabel],
            alpha=patch_alpha
        )
        plt.text(
            (r_width+l_width/2) * xMax,
            rightWidths[rightLabel]['bottom'] + 0.5*rightWidths[rightLabel]['right'],
            rightLabel,
            {'ha': 'center', 'va': 'center'},
            fontsize = fontsize
        )

    # Plot strips
    for leftLabel in left_labels:
        for rightLabel in right_labels:
            labelColor = leftLabel
            if len(data[(data.left == leftLabel) & (data.right == rightLabel)]) >0 :
                # Create array of y values for each strip, half at let value,
                # half at right
                ys_d = np.array(50 * [leftWidths[leftLabel]['bottom']] + 50 * [rightWidths[rightLabel]['bottom']])
                ys_d = np.convolve(ys_d, 0.05*np.ones(20), mode='valid')
                ys_d = np.convolve(ys_d, 0.05*np.ones(20), mode='valid')
                ys_u = np.array(50 * [leftWidths[leftLabel]['bottom'] + ns_l[leftLabel][rightLabel]] + 50 * [rightWidths[rightLabel]['bottom'] + ns_r[leftLabel][rightLabel]])
                ys_u = np.convolve(ys_u, 0.05 * np.ones(20), mode='valid')
                ys_u = np.convolve(ys_u, 0.05 * np.ones(20), mode='valid')

                # Update bottom edges at each label so next strip starts at the right place
                leftWidths[leftLabel]['bottom'] += ns_l[leftLabel][rightLabel]
                rightWidths[rightLabel]['bottom'] += ns_r[leftLabel][rightLabel]
                plt.fill_between(
                    np.linspace(0, xMax, len(ys_d)), ys_d, ys_u, alpha=link_alpha,
                    color = color_dict[labelColor]
                )
    plt.gca().axis('off')
    plt.gcf().set_size_inches(6, 6)
    if save != None:
        plt.savefig(save, bbox_inches='tight', dpi=150)
    if close_plot:
        plt.close()