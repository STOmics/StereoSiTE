#!/usr/bin/env python3
# _*_ coding: utf-8 _*_

# @Author: LiuXing liuxing2@genomics.cn 
# @Date: 2022-06-15 17:27:47 
# @Last Modified by:   LiuXing 
# @Last Modified time: 2022-06-15 17:27:47 

from email.policy import default
import sys, os

import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import cell2location
import scvi
from scvi import REGISTRY_KEYS
import matplotlib
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 #enables correct plotting of text
import seaborn as sns
from optparse import OptionParser

class Cell2location():
    """
    Run cell type deconvolution or annotation by cell2location, including reference and ST(Spatial Transcriptomics) data procession.

    Parameters
    ----------
    ref_file
        reference file path, which can be csv or anndata format.
    adata_vis_file
        ST file path, which will be deconvoluted or annotated
    out_dir
        output directory path that will be used to save result. Three subdirectory will be created:
        reference_signatures: save processed single-cell reference in csv format if the ref_file was given in anndata format
        cell2location_map: save deconvolution model and result in anndata format
        figures: save figures generated in the deconvolution or annotation process.
    bin_size
        the bin size of ST data in h5ad file, default = 50. If the given data was in single-cell resolution, please specify bin_size=1.
        the N_cells_per_location will be caculated based on this parameter. N_cell_per_location = (bin_size*500/10000)^2 if bins_size > 1 else 1
    gpu
        Load model on default GPU if available (if None or True), or index of GPU to use (if int), or name of GPU (if str), or use CPU (if False).
    """
    def __init__(self, ref_file: str, 
                 adata_vis_file: str,
                 out_dir: str = os.getcwd(),
                 bin_size: int = 50, 
                 gpu: int = 1):
        self.ref_file = ref_file
        self.adata_vis_file = adata_vis_file
        self.results_folder = out_dir
        self.N_cells_per_location = int((bin_size*500/10000)**2) if bin_size > 1 else 1 #calculate the cell number per bin based on bin size
        self.gpu = f"cuda:{gpu}" if gpu else False
        self.ref_run_name = f'{self.results_folder}/reference_signatures'
        self.run_name = f'{self.results_folder}/cell2location_map'
        self.figures = f'{self.results_folder}/figures'
        os.makedirs(self.ref_run_name, exist_ok=True)
        os.makedirs(self.run_name, exist_ok=True)
        os.makedirs(self.figures, exist_ok=True)
        sc.settings.figdir = self.figures

    def run_deconvolution(self):
        """
        run both reference processing and deconvolution with default parameters
        """
        if self.ref_file.endswith(".h5ad"):
            inf_aver = self.process_ref()
        elif self.ref_file.endswith(".csv"):
            inf_aver = pd.read_csv(self.ref_file, index_col=0)
        adata_vis  = self.process_vis(inf_aver)
        return adata_vis

    def process_ref(self,
                    batch_key: str = 'sample',
                    labels_key: str = 'cell_type',
                    max_epochs: int = 1500) -> pd.DataFrame:
        """
        process single-cell sequence reference data

        Parameters
        ----------
        batch_key
            specify the key for obtaining batch infomation. For example, if the reference data collected from different sample, reaction or exprement batch,
            data source should be markered with batch key, and the program will revise the batch effect.
            default = sample
        labels_key
            key name for getting cell type annotation information, default = cell_type.
        max_epochs
            maximal epochs for model training, default=1500.
        
        Returns
        ----------
            infered average gene expression vector of every cell type
        """
        sc.settings.figdir = self.figures
        inf_aver_file = f"{self.ref_run_name}/inf_aver.csv"
        if (os.path.exists(inf_aver_file)):
            inf_aver = pd.read_csv(inf_aver_file, index_col=0)
            return inf_aver
        elif (os.path.exists(f"{self.ref_run_name}/model.pt")):
            adata_ref_file = f"{self.ref_run_name}/sc.h5ad"
            adata_ref = sc.read_h5ad(adata_ref_file)
            mod = cell2location.models.RegressionModel.load(self.ref_run_name, adata_ref, use_gpu=self.gpu)
        else:
            adata_ref = sc.read_h5ad(self.ref_file)
            adata_ref.obs_names_make_unique()
            sc.pp.filter_cells(adata_ref, min_genes = 0)
            sc.pp.filter_genes(adata_ref, min_counts = 10)
            adata_ref.var['mt'] = adata_ref.var_names.str.startswith(('mt-', 'MT-'))
            sc.pp.calculate_qc_metrics(adata_ref, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
            sc.pl.violin(adata_ref, ['n_genes', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, save='ref_violin.png')
            sc.pl.scatter(adata_ref, x='total_counts', y='n_genes', save = 'ref_scatter.png')

            from cell2location.utils.filtering import filter_genes
            selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
            plt.savefig(f"{self.figures}/ref_filter_genes.png")
            plt.clf()
            #filter the object
            adata_ref = adata_ref[:, selected].copy()
            #prepare anndata for the regression model
            cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                                                   # reaction / sample / batch
                                                   batch_key = batch_key,
                                                   #batch_key = 'Sample',
                                                   # cell type, covariate used fir contructing signatures
                                                   labels_key = labels_key,
                                                   # multiplicative technical effects (platform, 3' vs 5', donor effect)
                                                   #categorical_covariate_keys = ['Experiment']
                                                   #categorical_covariate_keys = ['Method']
                                                  )
            # create the regression model
            from cell2location.models import RegressionModel
            mod = RegressionModel(adata_ref)
            #vies anndata_setup as a sanity check
            mod.view_anndata_setup()
            mod.train(max_epochs=max_epochs, use_gpu=self.gpu)

            # In this section, we export the estimated cell abundance (summary of the posterior distribution).
            adata_ref = mod.export_posterior(
                adata_ref, sample_kwargs={'num_samples':1000, 'batch_size': 2500, 'use_gpu': self.gpu}
            )

            #Save model
            mod.save(f"{self.ref_run_name}", overwrite=True)

            #Save anndata object with results
            adata_file = f"{self.ref_run_name}/sc.h5ad"
            adata_ref.write(adata_file)
            mod.plot_history(20)
            plt.savefig(f"{self.figures}/ref_train_history.png")
            plt.clf()
            inf_aver = mod.samples[f"post_sample_means"]["per_cluster_mu_fg"].T
            if "detection_y_c" in list(mod.samples[f"post_sample_means"].keys()):
                inf_aver = inf_aver * mod.samples[f"post_sample_means"]["detection_y_c"].mean()
            aver = mod._compute_cluster_averages(key=REGISTRY_KEYS.LABELS_KEY)
            aver = aver[mod.factor_names_]
            plt.hist2d(np.log10(aver.values.flatten() + 1),
                        np.log10(inf_aver.flatten() + 1),
                        bins = 50,
                        norm = matplotlib.colors.LogNorm(),
            )
            plt.xlabel("Mean expression for every gene in every cluster")
            plt.ylabel("Estimated expression for every gene in every cluster")
            plt.savefig(f"{self.figures}/ref_train_QC.png")
            plt.clf()
            #export estimated expression in each cluster
        if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
            inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                                          for i in adata_ref.uns['mod']['factor_names']]].copy()
        else:
            inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                        for i in adata_ref.uns['mod']['factor_names']]].copy()
        inf_aver.columns = adata_ref.uns['mod']['factor_names']
        inf_aver.to_csv(inf_aver_file)
        return inf_aver

    def process_vis(self, 
                    inf_aver: pd.DataFrame, 
                    max_epochs: int = 5000,
                    batch_size: int = 90000,
                    anno: str = 'cell2loc_anno',
                    spot_size: int = 70) -> anndata:
        """
        deconvolute ST data based on the processed SC reference

        Parameters
        ----------
        inf_aver
            infered average gene expression vector of every cell type
        max_epochs
            maximal epochs for model training, default=5000.
        batch_size
            batch size that determines the amount of data loaded into the gpu memory, default=90000. 
            If the running report outOfMemory error, reduce the batch_size can help to resolve.
        anno
            label key that will be used to store annotation result, default=cell2loc_anno.
        spot_size
            specify spot size when draw bins or cells in space.
        Returns
        -----------
            anndata with deconvolution and annotation result.    
        """
        sc.settings.figdir = self.figures
        # find shared genes and subset both anndata and reference signatures
        adata_vis = anndata.read(self.adata_vis_file)
        if adata_vis.raw != None:
            adata_vis = adata_vis.raw.to_adata()
        adata_vis.raw = adata_vis
        #find mitochondria-encoded (MT) genes
        adata_vis.var['MT_gene'] = adata_vis.var_names.str.startswith(('mt-', 'MT-'))

        # remove MT genes for spatial mapping (keeping their counts in the object)
        adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()
        adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]
        intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
        adata_vis = adata_vis[:, intersect].copy()
        inf_aver = inf_aver.loc[intersect, :].copy()

        #prepare anndata for cell2location model
        cell2location.models.Cell2location.setup_anndata(adata=adata_vis)

        # create and train the model
        mod = cell2location.models.Cell2location(
            adata_vis, cell_state_df=inf_aver,
            # the expected average cell abundance: tissue-dependent
            # hyper-prior which can be estimated from paired histology:
            N_cells_per_location=self.N_cells_per_location,
            # hyperparameter controlling normalisation of
            # within-experiment variation in RNA detection:
            detection_alpha=20
        )
        mod.view_anndata_setup()
        
        if mod.adata.n_obs < batch_size:
            batch_size = None 
        mod.train(max_epochs=max_epochs,
                # train using full data (batch_size=None)
                batch_size=batch_size,
                # use all data points in training because
                # we need to estimate cell abundance at all locations
                train_size=1,
                use_gpu=self.gpu)

        # In this section, we export the estimated cell abundance (summary of the posterior distribution).
        adata_vis = mod.export_posterior(
            adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': 10000, 'use_gpu': self.gpu}
        )

        # Save model
        mod.save(f"{self.run_name}", overwrite=True)
        mod.plot_history(20)
        plt.legend(labels=['full data training'])
        plt.savefig(f"{self.figures}/vis_train_history.png")
        plt.clf()

        use_n_obs = 1000
        ind_x = np.random.choice(mod.adata_manager.adata.n_obs, np.min((use_n_obs, mod.adata.n_obs)), replace=False)
        mod.expected_nb_param = mod.module.model.compute_expected(
            mod.samples[f"post_sample_means"], mod.adata_manager, ind_x=ind_x
        )
        x_data = mod.adata_manager.get_from_registry(REGISTRY_KEYS.X_KEY)[ind_x, :]
        x_data = np.asarray(x_data.toarray())
        mod.plot_posterior_mu_vs_data(mod.expected_nb_param["mu"], x_data)
        plt.savefig(f"{self.figures}/vis_QC.png")
        plt.clf()

        # Save anndata object with results
        adata_file = f"{self.run_name}/sp.h5ad"
        #adata_vis.write(adata_file)
        # add 5% quantile, representing confident cell abundance, 'at least this amount is present',
        # to adata.obs with nice names for plotting
        adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
        cellList=list(set(inf_aver.columns) & set(adata_vis.obs.columns))
        with mpl.rc_context({'axes.facecolor':  'black', 'figure.figsize': [5, 5]}):
            sc.pl.spatial(adata_vis, img_key="hires", color=cellList, spot_size=spot_size, vmin=0, vmax='p99.2', 
                        cmap='magma', save="cell_abundance.png")

        adata_vis.obs[anno] = adata_vis.obs[cellList].idxmax(axis=1)
        
        # Save anndata object with results
        self._plotCells(anno, adata_vis, spot_size=spot_size)
        #self._plotMarkerGenes(anno, adata_vis)
        adata_vis.write(adata_file)
        return adata_vis

    def _plotCells(self,
                   anno: str,
                   adata_vis: anndata,
                   spot_size: int = 70):
        """
        draw the annotation result
        """
        import math
        from tqdm import tqdm
        rcParams['figure.figsize'] = 10, 10
        sc.pl.spatial(adata_vis, img_key="hires", color=anno, spot_size=spot_size, save=f"{anno}.png")

        cellsCount = adata_vis.obs[anno].value_counts()
        cellsRate = cellsCount/cellsCount.sum()*100
        cellsdf = pd.concat([cellsCount, cellsRate], axis=1)
        cellsdf.columns = ['count', 'rate']
        cellsdf.to_csv(f"{self.run_name}/{anno}_cell_count.tsv", sep="\t")
        plt.figure(figsize=(8, 8))
        cellsdf['count'].plot(kind = 'barh')
        i = 0
        sum = cellsdf['count'].sum()
        for _, v in cellsdf['rate'].items():
            plt.text(sum*(v+2)/100, i, '%.2f' % v, ha='center', va='bottom', fontsize=11)
            i+=1
        plt.savefig(f"{self.figures}/{anno}_cell_count.png")
        plt.clf()
        plotcol = 4
        cells = adata_vis.obs[anno].unique()
        plotrow = math.ceil(len(cells)/plotcol)
        figSize = (16, plotrow*4)
        fig = plt.figure(figsize=figSize,dpi=100)
        for j in tqdm(range(len(cells))):
            cell = cells[j]
            i = j+1
            row = int(i/plotcol)
            col = i - row*plotcol
            ax = plt.subplot(plotrow, plotcol, i)
            sc.pl.spatial(adata_vis, img_key="hires", color=anno, groups = [cell], spot_size=spot_size, show=False, ax = ax, title ="{0} ({1})".format(cell, cellsCount.loc[cell]), legend_loc=None)
            ax.set_xlabel("")
            ax.set_ylabel("")
        plt.savefig(f"{self.figures}/{anno}_split.png")
        plt.clf()

def main():
    """
    This program can be used to process scRNAseq reference and annotate ST(spatial transcriptomics) data
    %prog [options]
    """
    parser = OptionParser(main.__doc__)
    parser.add_option("-r", "--reference", action = "store", type = "str", dest = "reference", help = "reference file path, can be h5ad or csv format file.")
    parser.add_option("-i", "--vis", action = "store", type = "str", dest = "vis", help = "input stereo-seq data in h5ad format.")
    parser.add_option("-o", "--outDir", action = "store", type = "str", dest = "outDir", help = "output directory path.")
    parser.add_option("-g", "--gpu", action = "store", type = "int", default = 1, dest = "gpu", help = "give cuda gpu name that will be used to run this program. default=1")
    parser.add_option("--bin_size", action="store", type = int, default = 50, dest = "bin_size", help = "bin size of the given ST data. default=50")

    opts, args = parser.parse_args()

    if (opts.reference == None or opts.vis == None or opts.outDir == None):
        sys.exit(not parser.print_help())

    cell2location = Cell2location(opts.reference, opts.vis, opts.outDir, bin_size=opts.bin_size, gpu = opts.gpu)
    cell2location.run_deconvolution()

if __name__ == "__main__":
    main()
