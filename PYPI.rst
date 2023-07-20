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

.. _scanpy: https://scanpy.readthedocs.io/en/stable/
.. _anndata: https://anndata.readthedocs.io/en/stable/
.. _squidpy: https://squidpy.readthedocs.io/en/stable/
.. _stereosite: https://github.com/STOmics/stereosite 