#!/usr/bin/env python
# coding: utf-8

# Description
# A function to create a joint adata object from (1) scRNA-seq adata, and (2) scATAC-seq adata objects.
# adata.X : cells x (genes-RNA + genes-ATAC)
# adata.obs : cell metadata
# adata.var : gene metadata (modality, for example)
# adata.obsm : ["X_umap_RNA"], [X_umap_ATAC], respectively
# Notes.
# 1) We assume that the two objects have the exact same group of cells
# 2) the order of the cells don't have to be the same. We will sort the cell orders.
# 3) the counts are log-normalized per each modality, thus should be taken with care, only for visualization in here.
# 4) Currently, exCellxgene breaks when it tries to compute differential gene expression, etc. This is on Alec's hands.
# 5) We assume that the datasets follow the standardization used for Neurips 2021 Open Problems dataset. 
# This should be explicitly clarified later on, regarding the field names, etc.


# import packages
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scanpy import AnnData
from scipy import sparse


#def join_multimodal_adatas_RNA_ATAC(
#    RNA_path: Union[str, PathLike, Iterator[str]],
#    ATAC_path: Union[str, PathLike, Iterator[str]],
#    output_path: Union[str, PathLike, Iterator[str]],) -> AnnData:
def join_multimodal_adatas_RNA_ATAC(RNA_path, ATAC_path, output_path):

    # import the datasets (RNA and ATAC)
    adata_GEX = sc.read_h5ad(RNA_path)
    adata_ATAC = sc.read_h5ad(ATAC_path)

    # Switch the peaks back to gene activity
    adata_ATAC_genes = ad.AnnData(X = adata_ATAC.obsm['gene_activity'],
                                 obs = adata_ATAC.obs, 
                                 var = pd.DataFrame(index=adata_ATAC.uns['gene_activity_var_names']))
    adata_ATAC_genes.obsm["X_umap"] = adata_ATAC.obsm["umap"]
    adata_ATAC_genes

    # Add the assay name on the adata.var_names
    # change the var names in adata_ATAC_genes
    adata_ATAC_genes.var_names = adata_ATAC_genes.var_names + "-ATAC"
    adata_ATAC_genes.var_names

    # change the var names in adata_GEX
    adata_GEX.var_names = adata_GEX.var_names + "-RNA"
    adata_GEX.var_names 

    # Step 0. Sort the indices for both adata objects (in case the cell_id orders are different between adata_GEX and adata_ATAC_genes)
    # adata_GEX

    # 1) adata.X
    # log-normalize the raw counts in adata.X
    # recover the raw counts
    adata_GEX.X = adata_GEX.layers["counts"].copy()
    adata_GEX.X.todense()
    # log-normalize the counts for RNA object
    sc.pp.normalize_total(adata_GEX, target_sum=1e4)
    sc.pp.log1p(adata_GEX)

    # create a dataframe with obs and var as indicex and columns
    df_counts_RNA = pd.DataFrame(adata_GEX.X.todense(), 
                                index=adata_GEX.obs_names,
                                columns=adata_GEX.var_names)

    # sort the index
    df_counts_RNA = df_counts_RNA.sort_index()
    adata_GEX.X = sparse.csr_matrix(df_counts_RNA.to_numpy())

    # 2) adata.obsm
    df_umap = pd.DataFrame(adata_GEX.obsm["X_umap"],
                    index = adata_GEX.obs_names)
    df_umap = df_umap.sort_index()
    adata_GEX.obsm["X_umap"] = df_umap.to_numpy()

    # 3) adata.obs
    adata_GEX.obs = adata_GEX.obs.sort_index()
    #sc.pl.umap(adata_GEX, color="cell_type")

    # adata_ATAC_genes

    # 1) adata.X
    # log-normalize the raw counts in adata.X -> ATAC data is not integer counts (or we never saved it separately)
    # recover the raw counts
    # adata_ATAC_genes.X = adata_ATAC_genes.layers["counts"].copy()
    # adata_ATAC_genes.X.todense()
    # log-normalize the counts for RNA object
    # sc.pp.normalize_total(adata_GEX, target_sum=1e4)
    # sc.pp.log1p(adata_GEX)

    # create a dataframe with obs and var as indicex and columns
    df_counts_ATAC = pd.DataFrame(adata_ATAC_genes.X.todense(), 
                                index=adata_ATAC_genes.obs_names,
                                columns=adata_ATAC_genes.var_names)

    # sort the index
    df_counts_ATAC = df_counts_ATAC.sort_index()
    adata_ATAC_genes.X = sparse.csr_matrix(df_counts_ATAC.to_numpy())

    # 2) adata.obsm
    df_umap = pd.DataFrame(adata_ATAC_genes.obsm["X_umap"],
                    index = adata_ATAC_genes.obs_names)

    df_umap = df_umap.sort_index()
    adata_ATAC_genes.obsm["X_umap"] = df_umap.to_numpy()

    # 3) adata.obs
    adata_ATAC_genes.obs = adata_ATAC_genes.obs.sort_index()
    #sc.pl.umap(adata_ATAC_genes, color="cell_type")

    # Step1. Concatenate the two anndata objects
    # Alternatively, we can concatenate two anndata objects using ad.concat
    # Put the assay name for adata.var just in case we need to differentiate where the marker comes from (usually it shoudl be obvious!)
    adata_GEX.var['assay'] = 'RNA'
    adata_ATAC_genes.var['assay'] = 'ATAC'

    # NOTE: the order of the adata.obs should matter in terms of concatenation (which is weird), 
    # so we'd need to sort the orders of both adata objects before the concatenation (done in the previous step)

    # concatenate the two adata objects
    adata_joint = ad.concat([adata_GEX, adata_ATAC_genes], axis=1, merge="first")

    # transfer annotations from both modalities
    adata_joint.obs["cell_type_RNA"] = adata_ATAC_genes.obs["rna_ann"]
    adata_joint.obs["cell_type_ATAC"] = adata_GEX.obs["atac_ann"]

    # (optional)
    # remove unnecessary fields (or fields partially filled with NaNs, that cause errors in exCellxgene)
    # Tip. generalize this to remove any fields filled with NaNs.
    del adata_joint.obs["pseudotime_order_ATAC"]
    del adata_joint.obs["pseudotime_order_GEX"]
    # del adata_joint.obs["atac_ann"]
    # del adata_joint.obs["rna_ann"]
    del adata_joint.obsm["X_umap"]

    # save the UMAP embeddings from both RNA and ATAC
    adata_joint.obsm["X_umap_RNA"] = adata_GEX.obsm["X_umap"]
    adata_joint.obsm["X_umap_ATAC"] = adata_ATAC_genes.obsm["X_umap"]

    # save the anndata to the output_path (path name should include the h5ad filename as well)
    adata_joint.write_h5ad(output_path)

    return adata_joint