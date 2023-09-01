|PyPI| |Downloads| |stars| |Discourse| |Zulip|

StereoSiTE - Spatial Transcriptome Analysis in Python
======================================================

**StereoSiTE** is a package for the analysis and visualization of spatial transcriptome data.
It builds on top of `anndata`_, `scanpy`_ and `squidpy`_, from which it inherits modularity and scalability.
It provides analysis tools to dissect cellular neighborhood based on cell composition and quantitatively define cell-cell communication in spatial.

StereoSiTE's key applications
------------------------------

- Cellular Neighborhood (CN) clustering based on cell composition of each bin
- Spatial Cell Interaction Intensity (SCII) analysis

Citation
---------

If you use `stereosite`_ in your work, please cite the publication as follows:

    **StereoSiTE: A framework to spatially and quantitatively profile the cellular neighborhood organized iTME**

    Xing Liu, Chi Qu, Chuandong Liu, Na Zhu, Huaqiang Huang, Fei Teng, Caili Huang, Bingying Luo, Xuanzhu Liu, Yisong Xu, Min Xie, Feng Xi, Mei Li, Liang Wu, Yuxiang Li, Ao Chen, Xun Xu, Sha Liao, Jiajun Zhang

    bioRxiv 2022.12.31.522366; doi: https://doi.org/10.1101/2022.12.31.522366

Installation
-------------

Install StereoSiTE via PyPi by running:
::

    pip install stereosite

or via raw code by running:
::

    git clone https://github.com/STOmics/StereoSiTE.git

    cd StereoSiTE

    python setup.py install

Run examples
------------

**Transfer spatial gene expression matrix (gem) to anndata**
::
    from stereosite.read.gem import Gem_Reader
    gem_reader = Gem_Reader(gem_file)

    #without cell mask, bin_size must be specified. Here bin_size=200
    adata = gem_reader.gem2anndata(200) 

    #or with cell mask, gem will be transfered to anndata in single-cell resolution.
    adata = gem_reader.gem_with_cellmask_2anndata(mask_file)

**Cellular Neighborhood (CN)**
::
    #calculate CN with cell type deconvolution result of squared bin data
    from stereosite.cn.cellneighbor import cn_deconvolve
    import anndata
    adata = anndata.read(adata_file)
    cn_deconvolve(adata, use_rep='q05_cell_abundance_w_sf') # use_rep specify matrix used to calculate cell composition of every bin

    #or calculate CN with annotated cell bin data
    from stereosite.cn.cellneighbor import cn_cellbin
    import anndata
    adata = anndata.read(adata_file)
    cn_cellbin(adata, 400, n_neighbors = 20, resolution = 0.4, min_dist = 0.1)

    #CN result visualization
    from stereosite.plot.cellneighbor import umap, heatmap, spatial
    spatial(adata, spot_size=20)
    umap(adata)
    heatmap(adata)

**Spatial Cell Interaction Intensity (SCII)**
::
    #calculate SCII with cell type annotated data. Choose LR database based on the sample type. CellChatDB provide mouse and human database
    #separatly. CellphoneDB provide only human database, but also can be used to analyse mouse data by transfer mouse gene into homologous 
    #human gene, which will be automaticaly done by the software.
    from stereosite.scii import intensities_count
    interactiondb_file = "./datasets/LR_database/CellChatDB.mouse.csv"
    scii_dict = intensities_count(adata, interactiondb_file, distance_threshold = 50, anno = 'cell_type')

    #SCII result visualization
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

    #show intensity of interaction meidated by specific LR pair and between specific cells in situ
    from stereosite.plot.intensity import intensity_insitu
    cells = ['Non-immune cells', 'M2-like']
    genes = ['Ptprc', 'Mrc1']
    intensity_insitu(adata, cells, genes, radius = 50)

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








