import os, sys
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import datasets
from sklearn.preprocessing import minmax_scale
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
    if isinstance(data, pd.Series):
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
            v = cc_tme_link[lr, tme]
            if v <= 0:
                continue
            left.append(f"{model_names[0]} {lr}")
            mid_l.append(f"{model_names[2]} {tme}")
            left_weight.append(v)
            mid_l_weight.append(v)
        for cc in range(ccs):
            v = lr_tme_link[cc, tme]
            if v <= 0:
                continue
            right.append(f"{model_names[1]} {cc}")
            mid_r.append(f"{model_names[2]} {tme}")
            right_weight.append(v)
            mid_r_weight.append(v)
    #mid_l_weight = tmes*ccs*[cc_tme_link.sum()/(tmes*ccs)]
    #mid_r_weight = tmes*lrs*[lr_tme_link.sum()/(tmes*lrs)]
    
    #dataFrame = pd.DataFrame({'left': left, 'right': right, 'mid_l': mid_l, 'mid_r': mid_r, 
    #                      'left_weight': left_weight, 'right_weight': right_weight, 
    #                          'mid_l_weight': mid_l_weight, 'mid_r_weight': mid_r_weight}, index=range(len(left)))
    data_left = pd.DataFrame({'left': left, 'mid_l': mid_l,  
                            'left_weight': left_weight, 
                              'mid_l_weight': mid_l_weight, 
                             }, index=range(len(left)))
    data_right = pd.DataFrame({'right': right, 'mid_r': mid_r, 
                          'right_weight': right_weight, 
                              'mid_r_weight': mid_r_weight
                             }, index=range(len(right)))
    return data_left, data_right

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
    dataFrame = pd.DataFrame({'left': left, 'right': right, 'left_weight': leftWeight, 'right_weight': rightWeight}, index=range(len(left)))
    return dataFrame

def lr_pathway_dict_generate(interactiondb_file:str,
                             lr_label:str='interaction_name',
                             pathway_label:str='pathway_name',
                            ):
    interactiondb =pd.read_csv(interactiondb_file)
    lr_pathway_dict = dict(zip(interactiondb[lr_label].values, interactiondb[pathway_label].values))
    return lr_pathway_dict

def interaction_matrix_decomposition(interaction_matrix, interactiondb_file,
                                     components:int=5,
                                     W_filter:float=0.5, H_filter:float=0.5,
                                     cell_level:int=0):
    filter_matrix = interaction_matrix.copy()
    if cell_level==1:
        filter_matrix.columns = filter_matrix.columns.swaplevel(i=0, j=1)
    cellchatDB = pd.read_csv(interactiondb_file)
    lr_pathway_dict = lr_pathway_dict_generate(interactiondb_file)
    sender_df = pd.DataFrame(index = [lr_pathway_dict[f'{x[0]}_{x[1]}'] if f'{x[0]}_{x[1]}' in lr_pathway_dict.keys() \
                                      else cellchatDB[(cellchatDB['ligand'] == x[0]) & (cellchatDB['receptor'] == x[1])]['pathway_name'].values[0] \
                                      for x in filter_matrix.index])
    cells = filter_matrix.columns.get_level_values(cell_level).unique()
    for cell in cells:
        sender_df[cell] = filter_matrix[cell].sum(axis=1).values
    sender_df = sender_df.groupby(sender_df.index).sum()
    from sklearn.decomposition import NMF
    model = NMF(n_components = components, init='random', random_state=0)
    W = model.fit_transform(sender_df.T.values)
    H = model.components_

    normalized_W = minmax_scale(W, axis=1)
    normalized_W[normalized_W < W_filter] = 0

    df_W = pd.DataFrame(normalized_W, index = sender_df.T.index)
    left = []
    right = []
    left_weight = []
    right_weight = []
    for cell, row in df_W.iterrows():
        for i, v in row.items():
            if v == 0:
                continue
            left.append(cell)
            right.append(f"Pattern {i+1}")
            left_weight.append(v)
            right_weight.append(v)
    sankey_df_W = pd.DataFrame({'left': left, 'right': right, 'left_weight': left_weight, 'right_weight': right_weight}, index=range(len(left)))

    normalized_H = minmax_scale(H, axis=0)
    normalized_H[normalized_H < H_filter] = 0
    df_H = pd.DataFrame(normalized_H, columns = sender_df.T.columns)
    left = []
    right = []
    left_weight = []
    right_weight = []
    average_v = normalized_H.sum()/normalized_H.shape[0]
    for i, row in df_H.iterrows():
        for pathway, v in row.items():
            if v == 0:
                continue
            left.append(f"Pattern {i+1}")
            right.append(pathway)
            left_weight.append(v)
            right_weight.append(v)
    sankey_df_H = pd.DataFrame({'left': left, 'right': right, 'left_weight': left_weight, 'right_weight': right_weight}, index=range(len(left)))
    return sankey_df_W, sankey_df_H

def sankey_3d(data_l:pd.DataFrame, data_r:pd.DataFrame, 
              cmap='tab20', left_labels = None, right_labels = None, mid_labels = None,
              aspect=3, fontsize=5, save=None, close_plot=False, patch_alpha:float=0.99, link_alpha:float=0.4,
              interval:float=0.005, module_color=(0.4980392156862745, 0.4980392156862745, 0.4980392156862745),
              ):
    '''
    Make Sankey Diagram showing flow:  ligand-receptor model <--- TME ---> cell-cell model
    
    Inputs:
        data_l: pandas.dataFrame.
            Contains columns left, mid_l, left_weight, mid_l_weight
        data_r: pandas.dataFrame.
            Contains columns right, mid_r, right_weight, mid_r_weight
        cmap: str|dict.
            Define colors of each patch. User can set matplotlib's colormap (e.g. viridis, jet, tab10) or label_name -> color dict (e.g. dict(A="red", B="blue", C="green", ...)).
        left_labels: list[str]|array[str].
            Order of the left labels in the diagram
        right_labels: list[str]|array[str].
            Order of the right labels in the diagram
        mid_labels: list[str]|array[str].
            Order of the middle labels in the diagram
        aspect: float.
            Vertical extent of the diagram in units of horizontal extent
        fontsize: float.
            Fontsize of patch label text
        save: str.
            If the figure file name was given, the sankey figure will be stored in it.
        interval: float.
            Distance between two adjacent patchs = interval * vertical length of all patchs.
        module_color: tuple[float, float, float].
            Define the color of left and right patchs representing modules. default=(0.4980392156862745, 0.4980392156862745, 0.4980392156862745).
        
    Output:
        None
    '''
    plt.figure(dpi=140)
    plt.rc('text', usetex=False)
    plt.rc('font', family='serif')

    if len(data_l[(data_l.left.isnull()) | (data_r.right.isnull()) | (data_l.mid_l.isnull()) | (data_r.mid_r.isnull())]):
        raise NullsInFrame('Sankey graph dose not support null values.')

    #Identify all labels that appear 'left' or 'right'
    allLabels = pd.Series(np.r_[data_r.mid_r.unique(), data_l.left.unique(), data_r.right.unique(), data_l.mid_l.unique()]).unique()

    #Identify labels
    if left_labels == None:
        left_labels = data_l.left.unique()[::-1]
    else:
        _check_data_matches_labels(left_labels, data_l['left'])
    if mid_labels == None:
        mid_labels = data_r.mid_r.unique()[::-1]
    else:
        _check_data_matches_labels(mid_labels, data_l['mid_l'])
    if right_labels == None:
        right_labels = data_r.right.unique()[::-1]
    else:
        _check_data_matches_labels(right_labels, data_r['right'])

    if isinstance(cmap, str):
        color_dict = {}
        colorPalette = sns.color_palette(cmap, len(mid_labels))
        for i, label in enumerate(mid_labels):
            color_dict[label] = colorPalette[i]
        for label in left_labels:
            color_dict[label] = module_color
        for label in right_labels:
            color_dict[label] = module_color
    elif isinstance(cmap, dict):
        color_dict=cmap
        missing = [label for label in allLabels if label not in color_dict.keys()]
        if missing:
            msg = "The cmap parameter is missing values for the following labels: {}".format(', '.join(missing))
            raise ValueError(msg)
    else:
        raise ValueError("cmap must be string representing the matplotlib's colormap or dict")

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
            midLDict[leftLabel] = data_l[(data_l.mid_l == midLabel) & (data_l.left == leftLabel)].mid_l_weight.sum()
            leftDict[leftLabel] = data_l[(data_l.mid_l == midLabel) & (data_l.left == leftLabel)].left_weight.sum()
        ns_m_l[midLabel] = midLDict
        ns_l[midLabel] = leftDict
        for rightLabel in right_labels:
            midRDict[rightLabel] = data_r[(data_r.mid_r == midLabel) & (data_r.right == rightLabel)].mid_r_weight.sum()
            rightDict[rightLabel] = data_r[(data_r.mid_r == midLabel) & (data_r.right == rightLabel)].right_weight.sum()
        ns_m_r[midLabel] = midRDict
        ns_r[midLabel] = rightDict

    midLWidths = defaultdict()
    for i, midLabel in enumerate(mid_labels):
        myD = {}
        myD['mid'] = data_l[data_l.mid_l == midLabel].mid_l_weight.sum()
        if i == 0:
            myD['bottom'] = 0
            myD['top'] = myD['mid']
        else:
            myD['bottom'] = midLWidths[mid_labels[i - 1]]['top'] + interval*data_l.mid_l_weight.sum()
            myD['top'] = myD['bottom'] + myD['mid']
            topEdge = myD['top']
        midLWidths[midLabel] = myD
    midRWidths = defaultdict()
    for i, midLabel in enumerate(mid_labels):
        myD = {}
        myD['mid'] = data_r[data_r.mid_r == midLabel].mid_r_weight.sum()
        if i == 0:
            myD['bottom'] = 0
            myD['top'] = myD['mid']
        else:
            myD['bottom'] = midRWidths[mid_labels[i - 1]]['top'] + interval*data_r.mid_r_weight.sum()
            myD['top'] = myD['bottom'] + myD['mid']
            topEdge = myD['top']
        midRWidths[midLabel] = myD

    # Determine positions of left label patches and total widths
    leftWidths = defaultdict()
    for i, leftLabel in enumerate(left_labels):
        myD = {}
        myD['left'] = data_l[data_l.left == leftLabel].left_weight.sum()
        if i == 0:
            myD['bottom'] = 0
            myD['top'] = myD['left']
        else:
            myD['bottom'] = leftWidths[left_labels[i - 1]]['top'] + interval*data_l.left_weight.sum()
            myD['top'] = myD['bottom'] + myD['left']
            topEdge = myD['top']
        leftWidths[leftLabel] = myD

    # Determine positions of right label patches and total widths
    rightWidths = defaultdict()
    for i, rightLabel in enumerate(right_labels):
        myD = {}
        myD['right'] = data_r[data_r.right == rightLabel].right_weight.sum()
        if i == 0:
            myD['bottom'] = 0
            myD['top'] = myD['right']
        else:
            myD['bottom'] = rightWidths[right_labels[i-1]]['top'] + interval*data_r.right_weight.sum()
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
            alpha=patch_alpha,
            edgecolor='k',
            linewidth=0.3,
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
            alpha=patch_alpha,
            edgecolor='k',
            linewidth=0.3,
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
            if len(data_l[(data_l.mid_l == midLabel)& (data_l.left == leftLabel)]) > 0:
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
            if len(data_r[(data_r.mid_r == midLabel) & (data_r.right == rightLabel)]) > 0 :
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
        plt.savefig(save, bbox_inches='tight', dpi=300)
    if close_plot:
        plt.close()

def sankey_2d(data:pd.DataFrame, cmap="tab20", left_labels:list=None, right_labels:list=None,
              aspect:float=4, fontsize=4, save=None, close_plot=False, patch_alpha:float=0.99, link_alpha:float=0.65,
              interval:float=0.0, strip_color='left', ax:mpl.axes.Axes=None, dpi:float=300
              ):
    '''
    Make Sankey Diagram
    
    Inputs:
        data: pandas.DataFrame.
            contains columns left, right, mid_l, mid_r, left_weight, right_weight, mid_l_weight, mid_r_weight
        cmap: str|dict.
            Define colors of each patch. User can set matplotlib's colormap (e.g. viridis, jet, tab10) or label_name -> color dict (e.g. dict(A="red", B="blue", C="green", ...)).
        left_labels: list[str] | array[str].
            order of the left labels in the diagram
        right_labels: list[str] | array[str].
            order of the right labels in the diagram
        aspect: float.
            vertical extent of the diagram in units of horizontal extent
        
    Output:
        None
    '''  
    if ax==None:
        fig = plt.figure(dpi=dpi)
        ax = fig.add_subplot(1, 1, 1)

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

    if isinstance(cmap, str):
        color_dict = {}
        colorPalette = sns.color_palette(cmap, len(allLabels))
        for i, label in enumerate(allLabels):
            color_dict[label] = colorPalette[i]
    elif isinstance(cmap, dict):
        color_dict = cmap
    else:
        raise Exception("cmap must be string representing the matplotlib's colormap or dict")

    #Determine widths of individual strips
    from collections import defaultdict
    ns_l = defaultdict()
    ns_r = defaultdict()
    for leftLabel in left_labels:
        leftDict = {}
        rightDict = {}
        for rightLabel in right_labels:
            leftDict[rightLabel] = data[(data.left == leftLabel) & (data.right == rightLabel)].left_weight.sum()
            rightDict[rightLabel] = data[(data.left == leftLabel) & (data.right == rightLabel)].right_weight.sum()
        ns_l[leftLabel] = leftDict
        ns_r[leftLabel] = rightDict

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
    l_width = -0.3
    r_width = 1.3
    xMax = topEdge/aspect

    #Draw vertical bars on left and right of each label's section & print label
    for leftLabel in left_labels:
        ax.fill_between(
            [l_width*xMax, 0],
            2 * [leftWidths[leftLabel]['bottom']],
            2 * [leftWidths[leftLabel]['bottom'] + leftWidths[leftLabel]['left']],
            color = color_dict[leftLabel],
            alpha=patch_alpha,
            edgecolor='k',
            linewidth=0.3,
        )
        ax.text(
            l_width/2 * xMax,
            leftWidths[leftLabel]['bottom'] + 0.5*leftWidths[leftLabel]['left'],
            leftLabel,
            {'ha': 'center', 'va': 'center'},
            fontsize=fontsize
        )
    for rightLabel in right_labels:
        ax.fill_between(
            [xMax, r_width*xMax], 2 * [rightWidths[rightLabel]['bottom']],
            2 * [rightWidths[rightLabel]['bottom'] + rightWidths[rightLabel]['right']],
            color = color_dict[rightLabel],
            alpha=patch_alpha,
            edgecolor='k',
            linewidth=0.3,
        )
        ax.text(
            (r_width+l_width/2) * xMax,
            rightWidths[rightLabel]['bottom'] + 0.5*rightWidths[rightLabel]['right'],
            rightLabel,
            {'ha': 'center', 'va': 'center'},
            fontsize = fontsize
        )

    # Plot strips
    for leftLabel in left_labels:
        for rightLabel in right_labels:
            if strip_color == 'left':
                labelColor = leftLabel
            else:
                labelColor = rightLabel
            if len(data[(data.left == leftLabel) & (data.right == rightLabel)]) > 0 :
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
                ax.fill_between(
                    np.linspace(0, xMax, len(ys_d)), ys_d, ys_u, alpha=link_alpha,
                    color = color_dict[labelColor]
                )
    ax.axis('off')
    if save != None:
        plt.savefig(save, bbox_inches='tight', dpi=dpi)
    if close_plot:
        plt.close()

def sankey_pathway_decomposition(W:pd.DataFrame, H:pd.DataFrame,
                                 left_labels:list=None, right_labels:list=None, mid_labels:list=None,
                                 aspect:float=4, patch_alpha:float=0.99, link_alpha:float=0.65,
                                 interval:float=0.0,
                                 figsize:tuple=(12, 6), fontsize:float=9,
                                 cmap='tab20',
                                 dpi:int=300,
                                 save:str=None,
                                 ):
    '''
    Input:
        W: pd.DataFrame.
            The DataFrame generated from W matrix of decomposition result with shape (K, R). K represents cells while R represents patterns.
        H: pd.DataFrame.
            The DataFrame generated from H matrix of decomposition result with shape (R, N). R represents 
        cmap: str|dict.
            Define colors of each patch. User can set matplotlib's colormap (e.g. viridis, jet, tab10) or label_name -> color dict (e.g. dict(A="red", B="blue", C="green", ...)).
    Output:
        None
    '''

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize, layout="constrained")
    sankey_2d(W, left_labels=left_labels, right_labels=mid_labels,
              aspect=aspect, patch_alpha=patch_alpha, link_alpha=link_alpha,
              interval=interval,
              cmap=cmap, fontsize=fontsize, strip_color='right', 
              dpi=dpi, ax=ax1)
    sankey_2d(H, left_labels=mid_labels, right_labels=right_labels,
              aspect=aspect, patch_alpha=patch_alpha, link_alpha=link_alpha,
              interval=interval,
              cmap=cmap, fontsize=fontsize, strip_color='left', 
              dpi=dpi, ax=ax2)
    if save != None:
        fig.savefig(save, bbox_inches='tight', dpi=dpi)