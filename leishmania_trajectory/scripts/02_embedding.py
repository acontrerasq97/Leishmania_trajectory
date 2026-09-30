#!/usr/bin/env python3

"""
Generate low-dimensional representations of the
preprocessed Leishmania scRNA-seq dataset.

This script generates:

    - PCA
    - KNN graph
    - UMAP
    - Diffusion map
    - FLE
"""

from pathlib import Path

import numpy as np
import scanpy as sc

# ============================================================
# 1. CONFIGURATION
# ============================================================

BASE_DIR = Path('./leishmania_trajectory')

DATA_DIR = BASE_DIR / "data"
PROCESSED_DIR = DATA_DIR / "processed"

INPUT_H5AD = PROCESSED_DIR / "ExprMatrix.preprocessed.h5ad"
OUTPUT_H5AD = PROCESSED_DIR / "ExprMatrix.embedded.h5ad"

# Parameters
N_PCS = 14
N_NEIGHBORS = 14
RANDOM_STATE = 42

# ================================================================
# 2. LOAD PREPROCESSED DATA
# ================================================================

adata = sc.read_h5ad(INPUT_H5AD)

# ================================================================
# 3. PCA
# ================================================================

# Use only HVGs
adata_hvg = adata[:, adata.var["highly_variable"]].copy()

sc.pp.scale(adata_hvg, max_value=10)

sc.tl.pca(
    adata_hvg, 
    n_comps=N_PCS, 
    svd_solver="arpack", 
    random_state=RANDOM_STATE
    )

# Transfer PCA coordinates to original AnnData
adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"].copy()

# Transfer PCA metadata
adata.uns["pca"] = adata_hvg.uns["pca"].copy()

adata.obsm['X_pca'].shape


# ================================================================
# 4. NEIGHBOR GRAPH
# ================================================================

sc.pp.neighbors(
    adata, 
    n_neighbors=N_NEIGHBORS, 
    n_pcs=N_PCS, 
    random_state=RANDOM_STATE
    )

adata.obsm['X_umap'] = adata.obs[['x', 'y']].to_numpy().copy()

#sc.pl.umap(adata, color=['day'])

# ================================================================
# 5. DIFFUSION MAP
# ================================================================



sc.tl.diffmap(adata, n_comps=N_PCS)

# ================================================================
# 6. CREATE EMBEDDING DATASET
# ================================================================

adata.write_h5ad(OUTPUT_H5AD, compression="gzip")
