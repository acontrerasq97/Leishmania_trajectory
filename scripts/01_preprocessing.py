"""
Preprocessing of Leishmania scRNA-seq data for trajectory inference.

"""
from pathlib import Path

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
OUTPUT_H5AD = PROCESSED_DIR / "ExprMatrix.preprocessed.h5ad"
OUTPUT_METADATA = PROCESSED_DIR / "cell_metadata.csv"
GENE_SET_SCORES_PATH = PROCESSED_DIR / 'gene_set_scores.csv'
CELL_SETS_PATH = BASE_DIR / 'data' / 'metadata' / 'cell_sets.gmt'

# Parameters
#N_HVG = 2000 # Only if the method Seurat v3 will be used

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
# Support raw counts
adata.layers["counts"] = adata.X.copy()


# ---------------------------------------
# Method: Seurat V3
# ---------------------------------------

#sc.pp.highly_variable_genes(
#    adata,
#    n_top_genes=N_HVG,
#    flavor="seurat_v3",
#    layer="counts"
#)

#adata.var['highly_variable'].sum()

# Save the result in a custom column
#adata.var["hvg_seurat_v3"] = adata.var["highly_variable"]

# ---------------------------------------------------
# Method: Seurat V2 DOI: 10.1016/J.CELL.2019.01.006
# ---------------------------------------------------

# Scale all cells to add up to 10,000
sc.pp.normalize_total(adata, target_sum=1e4)

# Apply log(x + 1) to stabilize the variance
sc.pp.log1p(adata)

# Method: mean.var.plot of Seurat v2 -> flavor='seurat'
sc.pp.highly_variable_genes(
    adata,
    flavor='seurat',           # Activates the classic method mean.var.plot (Seurat v2)
    n_bins=20,                 # Divides the genes into 20 bins
    min_disp=1.0,              # Threshold of the Z-scored Fano Factor (Dispersion) at 1.0
    max_disp=float('inf'),     # No maximum limit for dispersion
    min_mean=0.1,              # Minimum expression threshold to avoid switched-off genes 
    max_mean=float('inf')      # No maximum average expression limit
)

n_hvg = adata.var['highly_variable'].sum()
n_hvg
#del adata.var["highly_variable"] # Cleaning default column

sc.pl.highly_variable_genes(adata) #Plot HVG


# ============================================================
# 5. CREATE PROCESSED DATASET
# ============================================================

# Saving processed dataset
adata.write_h5ad(OUTPUT_H5AD, compression="gzip")

# Metadata separately
adata.obs.to_csv(OUTPUT_METADATA)

# ============================================================
# FINAL SUMMARY
# ============================================================

print(f"Cells          : {adata.n_obs}")
print(f"Genes          : {adata.n_vars}")
print(f"HVGs           : {n_hvg}")

print("\nCells per day:")
print(
    adata.obs["day"]
    .value_counts()
    .sort_index()
)

print("\nOutput:")
print(OUTPUT_H5AD)