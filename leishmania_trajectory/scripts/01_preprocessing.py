#!usr/bin/env python

"""
Preprocessing of Leishmania scRNA-seq data for trajectory inference.

"""
import os
from pathlib import Path

import ipywidgets as widgets
import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc
import wot
from matplotlib import pyplot as plt

# ================================================================
# 1. CONFIGURATION
# ================================================================ 

BASE_DIR = Path('./leishmania_trajectory')

DATA_DIR = BASE_DIR / 'data'/ 'principal'
PROCESSED_DIR = BASE_DIR / 'data'/ 'processed'
RESULTS_DIR = BASE_DIR / 'results'

#PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
#RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Input files
FULL_DS_PATH = BASE_DIR / 'data'/ 'principal' / 'ExprMatrix.h5ad'
UMAP_COORDS_PATH = BASE_DIR / 'data' / 'metadata' / 'umap_coords.txt'
CELL_DAYS_PATH = BASE_DIR / 'data' / 'metadata' /'cell_days.txt'
GENE_SETS_PATH = BASE_DIR / 'data' / 'metadata' /'gene_sets.gmx'

# Output files
VAR_DS_PATH = BASE_DIR / 'data' / 'metadata' / 'ExprMatrix.preprocessed.h5ad'
GENE_SET_SCORES_PATH = BASE_DIR / 'data' / 'metadata' / 'gene_set_scores.csv'
CELL_SETS_PATH = BASE_DIR / 'data' / 'metadata' / 'cell_sets.gmt'

# Parameters
N_HVG = 2000

TIME_ORDER = ["0", "2", "5", "9"]


# ================================================================
# 2. LOAD DATA
# ================================================================

adata = sc.read_h5ad(FULL_DS_PATH) # loading expression matrix
coord_df = pd.read_csv(UMAP_COORDS_PATH, index_col='id', sep='\t')
days_df = pd.read_csv(CELL_DAYS_PATH, index_col='id', sep='\t')

adata.obs = adata.obs.merge(coord_df, left_index=True, right_index=True)

unique_days = adata.obs['day'].unique()
unique_days = unique_days[np.isnan(unique_days) == False]

# plot visualization coordinates
fig, ax = plt.subplots(figsize=(10, 10))
ax.axis('off')

scatter = ax.scatter(
    adata.obs['x'], 
    adata.obs['y'], 
    c=adata.obs['day'],
    s=1, 
    marker='o', 
    edgecolors='none', 
    alpha=0.8, 
    cmap='viridis'
)

cb = fig.colorbar(scatter, ax=ax)
cb.ax.set_title('Day')

plt.tight_layout()
plt.show()
plt.close(fig)

# =========================================================
# 3. GENE-SET SCORES
# =========================================================

if GENE_SETS_PATH.exists():

    gs = wot.io.read_sets(GENE_SETS_PATH, adata.var.index.values)
    gene_set_scores_df = pd.DataFrame(index=adata.obs.index)

    for j in range(gs.shape[1]):
        gene_set_name = str(gs.var.index.values[j])
        result = wot.score_gene_sets(ds=adata, gs=gs[:, [j]], permutations=0, method='mean_z_score')
        gene_set_scores_df[gene_set_name] = result['score']

    gene_set_scores_df.to_csv(GENE_SET_SCORES_PATH, index_label='id')
else:

    print(
        "WARNING: gene_sets.gmx was not found. "
        "Skipping gene-set scoring."
    )

# ============================================================
# 4. HIGHLY VARIABLE GENES
# ============================================================


# Identifying highly variable genes
adata.layers["counts"] = adata.X.copy()

#Apply log(x + 1) to stabilize the variance
sc.pp.log1p(adata)
