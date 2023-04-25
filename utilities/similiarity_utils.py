#!/usr/bin/env python
# coding: utf-8

# Description
# A function to compute the similarity metrics between two groups of cells
# group1: cells by features 1 
# group2: cells by features 2

# Notes.
# We assume that the two objects have the exact same group of cells

# examples of the features: 
# 1) genes (RNA)
# 2) genes/peaks (ATAC)
# 3) PCs (from PCA)
# 4) LSIs (from ATAC)
# 5) list of markers (genes, peaks, PCs, etc.)

# import packages
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy import sparse
# import modules for computing the similarity metrics (between two cells)
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.metrics.pairwise import euclidean_distances
from scipy.stats import pearsonr

def compute_similarity_metrics(group1, group2, metric):
    if len(group1)!=len(group2):
        error("the lengths of group1 and group2 are different.")
    
    if metric=="cosine": # cosine similarity
        similarity = cosine_similarity(group1, group2).mean(axis=1)
    elseif metric=="eucilidean": # eucilidean distance
        similarity = euclidean_distances(group1, group2).mean(axis=1)
    elseif metric=="correlation":
        similiarity = pearsonr
    else:
        error("metric is not define. Make sure to check the input arguments.")

    return similarity
