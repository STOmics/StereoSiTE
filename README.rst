|PyPI| |Downloads| |stars| |Discourse| |Zulip|

StereoSiTE - Spatial Transcriptome Analysis in Python
======================================================

**StereoSiTE** is a package for the analysis and visualization of spatial transcriptome data.
It builds on top of `anndata`_, `scanpy`_ and `squidpy`_, from which it inherits modularity and scalability.
It provides analysis tools to dissect cellular neighborhood based on cell composition and quantitatively define cell-cell communication in spatial.

StereoSiTE's key applications
------------------------------

- Cellular Neighborhood (CN) clustering based on cellular composition of each bin
- Spatial Cell Interaction Intensity (SCII) analysis
- Tissue domain clustering based on intercellular interactions within each bin (SCIITensor)

Demo data
----------
`demo`_ (Fetch code: nFks)

Citation
---------

If you use `stereosite`_ in your work, please cite the publication as follows:

    **StereoSiTE: A framework to spatially and quantitatively profile the cellular neighborhood organized iTME**

    Xing Liu, Chi Qu, Chuandong Liu, Na Zhu, Huaqiang Huang, Fei Teng, Caili Huang, Bingying Luo, Xuanzhu Liu, Min Xie, Feng Xi, Mei Li, Liang Wu, Yuxiang Li, 
    Ao Chen, Xun Xu, Sha Liao, Jiajun Zhang, StereoSiTE: a framework to spatially and quantitatively profile the cellular neighborhood organized iTME, 
    GigaScience, Volume 13, 2024, giae078, https://doi.org/10.1093/gigascience/giae078

    **SCIITensor: A tensor decomposition based algorithm to construct actionable TME modules with spatially resolved intercellular communications**

    Huaqiang Huang, Chuandong Liu, Xin Liu, Jingyi Tian, Feng Xi, Mei Li, Guibo Li, Ao Chen, Xun Xu, Sha Liao, Jiajun Zhang, Xing Liu
    bioRxiv 2024.05.21.595103; doi: https://doi.org/10.1101/2024.05.21.595103


Installation
-------------

Install StereoSiTE via PyPi by running:

    >>> conda create -y -n stereosite python=3.9.12
    >>> conda activate stereosite
    >>> pip install stereosite

or via raw code by running:

    >>> git clone https://github.com/STOmics/StereoSiTE.git
    >>> cd StereoSiTE
    >>> python setup.py install

Run examples
-------------

**Transfer spatial gene expression matrix (gem|gef) to anndata**
::

    from stereosite.read.gem import Gem_Reader
    gem_file = "to/path/sn.gem.gz"
    gem_reader = Gem_Reader(gem_file)
    # without cell mask, bin_size must be specified. Here bin_size=200
    adata_bin200 = gem_reader.gem2anndata(200)
    adata_bin200_file = "to/path/sn_bin200.h5ad"
    adata_bin200.write(adata_bin200_file)
    # or with cell mask, gem will be transfered to anndata in single-cell resolution.
    mask_file = "to/path/sn_mask.tif"
    adata = gem_reader.gem_with_cellmask_2anndata(mask_file)

**Cellular Neighborhood (CN)**
::

    from stereosite.cn.deconvolution import Cell2location
    # Get the cellular composition of each square bin by deconvolution
    ref_file = "to/path/scCell_reference.csv"
    adata_file = "to/path/sn_bin200.h5ad"
    out_dir = "to/out/deconvolution/bin200"
    cell2loc = Cell2location(ref_file, adata_file, out_dir = out_dir, bin_size = 200, gpu = 0)
    cell2loc.run_deconvolution()
    # Analyze the Cellular Neighborhood based on deconvolution result
    from stereosite.cn.cellneighbor import cn_deconvolve
    import anndata
    adata_anno_file = "to/out/deconvolution/bin200/cell2location_map/sp.h5ad"
    adata = anndata.read(adata_anno_file)
    # use_rep specify matrix used to calculate cell composition of every bin
    cn_deconvolve(adata, use_rep='q05_cell_abundance_w_sf')
    # Or use the annotated cell bin data to dissect the CNs
    from stereosite.cn.cellneighbor import cn_cellbin
    import anndata
    adata = anndata.read(adata_anno_file)
    cn_cellbin(adata, 400, n_neighbors = 20, resolution = 0.4, min_dist = 0.1)
    # CN result visualization
    from stereosite.plot.cellneighbor import umap, heatmap, spatial
    spatial(adata, spot_size=20)
    umap(adata)
    heatmap(adata)

**Spatial Cell Interaction Intensity (SCII)**
::

    from stereosite.scii import intensities_count
    # The annotated cellbin or square bin at single-cell resolution data is required.
    # Choose LR database based on the sample type. CellChatDB provide mouse and human database separatly. 
    # CellphoneDB provide only human database, but also can be used to analyse mouse data by transfer mouse gene into homologous
    # human gene, which will be automaticaly done by the software.
    interactiondb_file = "./datasets/LR_database/CellChatDB.mouse.csv"
    scii_dict = intensities_count(adata, interactiondb_file,
                                     distance_threshold = 50, 
                                     anno = 'cell_type')
    # Or we can specify different distance_threshold for individual LR types
    scii_dict = intensities_count(adata, interactiondb_file,
                                     distance_threshold = {'Secreted Signaling': 200, 'ECM-Receptor': 200, 'Cell-Cell Contact': 30}, 
                                     anno = 'cell_type')
    # Or specify the distance_coefficient parameter to consider distance when caculating interaction intensity.
    # distance_coefficient=0 means distance would not influence the interaction intensity.
    scii_dict = intensities_count(adata, interactiondb_file,
                                     distance_threshold = {'Secreted Signaling': 200, 'ECM-Receptor': 200, 'Cell-Cell Contact': 30},
                                     distance_coefficient = {'Secreted Signaling': 1, 'ECM-Receptor': 0.1, 'Cell-Cell Contact': 0},
                                    anno = 'cell_type')
    # The interaction result can be writen into a pickle file, and can be re-loaded when you want to re-analyze it.
    import pickle
    os.makedirs("./out/scii", exist_ok=True)
    interaction_file = "./out/scii/interactions.pkl"
    with open(interaction_file, 'wb') as writer:
         pickle.dump(scii_dict, writer)
    with open(interaction_file, 'rb') as reader:
         scii_dict = pickle.load(reader)

    # SCII result visualization
    #filter scii result based on intensities value or pathway
    cell_pairs = [('cell1', 'cell2'), ('cell1', 'cell3'), ...] # cell pairs that will be remained
    gene_list = ['gene1', 'gene2', 'gene3', ...] # list of genes that will be filtered out
    filter_interaction = interaction_select(scii_dict, cell_pairs=cell_pairs, filter_genes=gene_list, intensities_range=(6000, 8000))
    pathway_names = ['EGF'] # Interactions of these pathways will be remained.
    pathway_interaction = interaction_pathway_select(scii_dict, pathway_name=pathway_names, interactiondb_file=interactiondb_file)
    # Users can combine the intensities filter with the pathway selection
    pathway_interaction = interaction_pathway_select(filter_interaction, pathway_name=pathway_names, interactiondb_file=interactiondb_file)

    # Visualize the interaction result by bubble plot
    from stereosite.plot.scii import ligrec
    import numpy as np
    ligrec(scii_dict,
         intensities_range=(50, np.inf),
         pvalue_threshold=0.05,
         alpha=1e-4,
         swap_axes=False,
         source_groups=["Non-immune cells", "M2-like", 'DC', 'Teff'],
         target_groups = ["M1-like", "M2-like", "Monocytes", "Teff", "CD8+ Tcells"],
         title=" ",
    )
    # Or visualize the selected interactions
    ligrec(pathway_interaction,
         pvalue_threshold=0.05,
         alpha=1e-4,
         swap_axes=False,
    )
    # Show spatial distribution of interaction intensity between specific cell pair meidated by specific LR pair.
    from stereosite.plot.intensity import intensity_insitu
    cells = ['Non-immune cells', 'M2-like']
    genes = ['Ptprc', 'Mrc1']
    intensity_insitu(adata, cells, genes, radius = 50, distance_coefficient=0.01, spot_size=5)

    # Visualize the interaction result by circle plot and graph
    from stereosite.plot import scii_circos, scii_net
    anno='cell_type'
    cell_colors = dict(zip(adata.obs[anno].cat.categories, adata.uns[f'{anno}_colors'])) # Define the color of sectors representing cells
    filter_matrix = filter_interaction['intensities'].fillna(0)
    scii_circos.cells_lr_circos(filter_matrix, cell_colors=cell_colors)
    pathway_matrix = pathway_interaction['intensities'].fillna(0)
    scii_circos.cells_lr_circos(pathway_matrix, cell_colors=cell_colors)
    scii_circos.cells_circos(filter_matrix)
    scii_circos.cells_circos(pathway_matrix)

    #Draw the network diagram based on the Graph generated previously.
    g1 = scii_net.lr_link_graph_generate(filter_matrix, cell_colors=cell_colors, reducer=6)
    scii_net.cell_lr_grap_plot(g1, figsize=10, vertex_label_size=6)
    g2 = scii_net.cell_graph_generate(filter_matrix, reducer=30, cell_colors=cell_colors)
    scii_net.cell_graph_plot(g2, vertex_label_size=8, figsize=5, edge_width=[0.5, 3])
    g3 = scii_net.lr_link_graph_generate(pathway_matrix, cell_colors=cell_colors, reducer=6)
    scii_net.cell_lr_grap_plot(g3, figsize=8, edge_width=[0.5, 2])
    g4 = scii_net.cell_graph_generate(pathway_matrix, reducer=15, cell_colors=cell_colors)
    scii_net.cell_graph_plot(g4, vertex_label_size=8, figsize=5, edge_width=[0.5, 3])
    

**SCIITensor -- single sample analysis**
::

    from stereosite import scii_tensor
    import anndata
    import pandas as pd
    import seaborn as sns
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import pickle
    import numpy as np
    import scanpy as sc
    # Generate interactiontensor object and evaluate the optimal combination of ranks
    adata = anndata.read(adata_anno_file)
    interactionDB = "./datasets/LR_database/CellChatDB.mouse.csv"
    sct = scii_tensor.InteractionTensor(adata, interactionDB=interactionDB)
    radius = {'Secreted Signaling': 100, 'ECM-Receptor': 100, 'Cell-Cell Contact': 30}
    scii_tensor.build_SCII(sct, radius=radius, window_size=200, anno_col='cell2loc_anno')
    scii_tensor.process_SCII(sct, zero_remove=True, log_data=True)
    reconstruction_errors = scii_tensor.evaluate_ranks(sct, use_gpu=True, device='cuda:1')
    # Visualize the reconstruction errors using line plot
    from stereosite.plot.scii_tensor import reconstruction_error_line
    reconstruction_error_line(reconstruction_errors, figsize=(4, 4))
    # Decompose the interaction tensor with optimal combination of ranks
    scii_tensor.SCII_Tensor(sct, rank=[15, 15, 15], device='cuda:0')
    with open("out/scii_tensor_res.pkl", "wb") as f:
        pickle.dump(sct, f)
    # spatial distribution of each TME module
    import scanpy as sc
    sc.pl.spatial(sct.adata, color='TME_module', img_key=None, spot_size=20)

    ## visualization of core matrix
    from stereosite.plot import sankey, scii_circos, scii_net
    #normalize the core matrix
    norm_core = scii_tensor.core_normalization(sct.core, feature_range=(0, 100))
    #process core matrix to generate dataFrame that will be used to draw sankey plot.
    left_df, right_df = sankey.core_process(norm_core)
    sankey.sankey_3d(left_df, right_df, link_alpha=0.5, interval=0.005)

    from stereosite.plot.scii_tensor import tme_core_heatmap, core_heatmap
    core_heatmap(norm_core) # 3D heatmap plot showing the core matrix
    tme_core_heatmap(sct.core, tme_number=1, figsize=(4, 4)) # 2D heatmap plot showing the result of one TME
    ## cell-cell factor heatmap
    from stereosite.scii_tensor import top_pair
    import seaborn as sns
    import matplotlib.pyplot as plt
    top_cc_pair = top_pair(sct, pair='cc', top_n=20)
    fig = sns.clustermap(top_cc_pair.T, cmap="Purples", standard_scale=0, metric='euclidean', method='ward', 
                          row_cluster=False, dendrogram_ratio=0.05, cbar_pos=(1.02, 0.6, 0.01, 0.3),
                          figsize=(4, 6),
                          )
    ## ligand-receptor factor heatmap
    top_lr_pair = top_pair(sct, pair='lr', top_n=20)
    fig = sns.clustermap(top_lr_pair.T, cmap="Purples", standard_scale=0, metric='euclidean', method='ward', 
                          row_cluster=False, dendrogram_ratio=0.05, cbar_pos=(1.02, 0.6, 0.01, 0.3),
                          figsize=(4, 6),
                          )

    ## visualize selected interactions using heatmap
    from stereosite.plot.scii_tensor import interaction_heatmap
    interactions = scii_tensor.interaction_select(sct,
                                                   tme_module=1,
                                                   cellpair_module=1,
                                                   lrpair_module=11, n_lr=15, n_cc=15)
    interaction_heatmap(interactions, figsize=(5, 3), vmax=50)
    ## visualize selected interactions using circle plot
    from stereosite.plot.scii_circos import cells_lr_circos, cells_circos, lr_circos
    cells = adata.obs['cell2loc_anno'].unique()
    cell_colors = dict(zip(adata.obs[anno].cat.categories, adata.uns[f'{anno}_colors'])) # Define the color of sectors representing cells
    scii_circos.cells_lr_circos(interaction_matrix, cells=cells, cell_colors=cell_colors, scii_tensor=True)
    #Draw the circos which only contains cell types and the links between them.
    scii_circos.cells_circos(interaction_matrix, cells, cell_colors=cell_colors, label_orientation='vertical', scii_tensor=True)
    #Draw circos which only contains ligand-receptor genes
    scii_circos.lr_circos(interaction_matrix, cells=cells, scii_tensor=True)

    #Draw the network diagram based on the Graph generated previously.
    from stereosite.plot.scii_net import lr_link_graph_generate
    g1 = lr_link_graph_generate(interaction_matrix, cells = cells, separator="_", cell_colors=cell_colors, scii_tensor=True)
    scii_net.cell_lr_grap_plot(g1, figsize=10,

**SCIITensor -- multiple sample analysis**
::

    adata_1 = anndata.read(adata_anno_file_1)
    ## decompose another sample data
    ## evaluate the optimal combination of ranks
    interactionDB = "./datasets/LR_database/CellChatDB.mouse.csv"
    sct_1 = scii_tensor.InteractionTensor(adata_1, interactionDB=interactionDB)
    radius = {'Secreted Signaling': 100, 'ECM-Receptor': 100, 'Cell-Cell Contact': 30}
    scii_tensor.build_SCII(sct_1, radius=radius, window_size=200, anno_col='cell2loc_anno')
    scii_tensor.process_SCII(sct_1, zero_remove=True, log_data=True)
    reconstruction_errors = scii_tensor.evaluate_ranks(sct_1, use_gpu=True, device='cuda:1')
    ## visualize the reconstruction errors using line plot
    from stereosite.plot.scii_tensor import reconstruction_error_line
    reconstruction_error_line(reconstruction_errors, figsize=(4, 4))
    scii_tensor.SCII_Tensor(sct_1, rank=(20, 20, 13), device='cuda:0')
    with open("out/scii_tensor_res_1.pkl", "wb") as f:
    pickle.dump(sct_1, f)
    ## merge decomposed matrices
    sct_merge = scii_tensor.merge_data([sct, sct_1], patient_id=['p1' ,'p2'])
    ## visualize the reconstruction errors
    from stereosite.plot.scii_tensor import reconstruction_error_line
    reconstruction_error_line(reconstruction_errors, figsize=(4, 4))
    scii_tensor.SCII_Tensor_multiple(sct_merge, rank=[15,15,10], device='cuda:1')    
    ## spatial distribution of meta-module
    sc.pl.spatial(sct_merge.adata[0], color=['TME_module', 'TME_meta_module'], img_key=None, spot_size=20)
    sc.pl.spatial(sct_merge.adata[1], color=['TME_module', 'TME_meta_module'], img_key=None, spot_size=20)
    #normalize the core matrix
    norm_core = scii_tensor.core_normalization(sct.core, feature_range=(0, 100))
    #process core matrix to generate dataFrame that will be used to draw sankey plot.
    left_df, right_df = sankey.core_process(norm_core)
    sankey.sankey_3d(left_df, right_df, link_alpha=0.5, interval=0.005)
    from stereosite.plot.scii_tensor import tme_core_heatmap, core_heatmap
    core_heatmap(norm_core) # 3D heatmap plot showing the core matrix
    tme_core_heatmap(sct.core, tme_number=1, figsize=(4, 4)) # 2D heatmap plot showing the result of one TME
    
    ## visualize selected interactions using heatmap
    from stereosite.plot.scii_tensor import interaction_heatmap
    interactions = scii_tensor.interaction_select_multiple(sct_merge,
                                                        tme_module=0, sample='p2', 
                                                        cellpair_module=0, 
                                                        lrpair_module=1, n_lr=15, n_cc=15)
    interaction_heatmap(interactions, figsize=(5, 3), vmax=10)
    ##visualize selected interactions using circle plot
    from stereosite.plot.scii_circos import cells_lr_circos, cells_circos, lr_circos
    cells = adata.obs['cell2loc_anno'].unique()
    cell_colors = dict(zip(adata.obs[anno].cat.categories, adata.uns[f'{anno}_colors'])) # Define the color of sectors representing cells
    scii_circos.cells_lr_circos(interaction_matrix, cells=cells, cell_colors=cell_colors, scii_tensor=True)
    #Draw the circos which only contains cell types and the links between them.
    scii_circos.cells_circos(interaction_matrix, cells, cell_colors=cell_colors, label_orientation='vertical', scii_tensor=True)
    #Draw circos which only contains ligand-receptor genes
    scii_circos.lr_circos(interaction_matrix, cells=cells, scii_tensor=True)

.. |stars| image:: https://img.shields.io/github/stars/STOmics/StereoSiTE?logo=GitHub&color=yellow
    :target: https://github.com/STOmics/StereoSiTE/stargazers

.. |PyPI| image:: https://img.shields.io/pypi/v/stereosite.svg
    :target: https://pypi.org/project/stereosite/
    :alt: PyPI

.. |Downloads| image:: https://static.pepy.tech/badge/stereosite
    :target: https://pepy.tech/project/stereosite
    :alt: Downloads

.. |Discourse| image:: https://img.shields.io/discourse/posts?color=yellow&logo=discourse&server=https%3A%2F%2Fdiscourse.scverse.org
    :target: https://discourse.scverse.org/
    :alt: Discourse

.. |Zulip| image:: https://img.shields.io/badge/zulip-join_chat-%2367b08f.svg
    :target: https://scverse.zulipchat.com
    :alt: Zulip

.. _scanpy: https://scanpy.readthedocs.io/en/stable/
.. _anndata: https://anndata.readthedocs.io/en/stable/
.. _squidpy: https://squidpy.readthedocs.io/en/stable/
.. _stereosite: https://github.com/STOmics/stereosite 
.. _demo: https://bgipan.genomics.cn/#/link/ilOA8JTgy7jKrNX4ZrOc








