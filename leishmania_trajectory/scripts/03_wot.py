#!/usr/bin/env python3

"""
Waddington-OT trajectory inference for the Leishmania
scRNA-seq dataset.

This script performs trajectory inference.
Biological validation is performed separately.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import wot

# ================================================================
# 1. CONFIGURATION
# ================================================================

BASE_DIR = Path("./leishmania_trajectory")

DATA_DIR = BASE_DIR / "data"
PROCESSED_DIR = DATA_DIR / "processed"

RESULTS_DIR = BASE_DIR / "results"
WOT_DIR = RESULTS_DIR / "wot"

WOT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

# ------------------------------------------------
# Input
# ------------------------------------------------

INPUT_H5AD = (
    PROCESSED_DIR /
    "ExprMatrix.preprocessed.h5ad"
)

GENE_SET_SCORES_PATH = (
    PROCESSED_DIR /
    "gene_set_scores.csv"
)


# ------------------------------------------------
# Output
# ------------------------------------------------

TMAP_PREFIX = WOT_DIR / "tmap_leishmania"

# ------------------------------------------------
# Wot Parameters
# ------------------------------------------------

EPSILON = 0.05

LAMBDA1 = 3

LAMBDA2 = 50

GROWTH_ITERS = 3

# WOT local PCA used internally to construct the
# expression-space cost matrix
LOCAL_PCA = 30

RANDOM_SEED = 42

# ------------------------------------------------
# Growth Model Parameters
# ------------------------------------------------

BETA_MAX = 1.7
BETA_MIN = 0.3

BETA_CENTER = 0.25
BETA_WIDTH = 0.5

DELTA_MAX = 1.7
DELTA_MIN = 0.3

DELTA_CENTER = 0.1
DELTA_WIDTH = 0.2

# ================================================================
# 2. LOAD DATA
# ================================================================

adata = sc.read_h5ad(INPUT_H5AD)

print(f"Cells : {adata.n_obs}")
print(f"Genes : {adata.n_vars}")


# ================================================================
# 3. LOAD GENE-SET SCORES
# ================================================================

gene_set_scores = pd.read_csv(GENE_SET_SCORES_PATH, index_col="id")

print(
    f"Gene-set score matrix: "
    f"{gene_set_scores.shape}"
)

# Check that all cells are present
missing_cells = adata.obs_names.difference(gene_set_scores.index)

extra_cells = gene_set_scores.index.difference(adata.obs_names)

if len(missing_cells) > 0:
    raise ValueError(
        f"Missing gene-set scores for "
        f"{len(missing_cells)} cells."
    )

if len(extra_cells) > 0:
    print(
        f"Warning: gene_set_scores contains "
        f"{len(extra_cells)} extra cells."
    )

# Check duplicated cell IDs
if gene_set_scores.index.duplicated().any():
    duplicated = gene_set_scores.index[
        gene_set_scores.index.duplicated()
    ]

    raise ValueError(
        f"Duplicated cell IDs found: "
        f"{duplicated.tolist()}"
    )

# Reorder
gene_set_scores = gene_set_scores.loc[adata.obs_names]

# Add scores to AnnData
for column in gene_set_scores.columns:

    adata.obs[column] = gene_set_scores[column].values

# ================================================================
# 8. GROWTH FUNCTIONS
# ================================================================

def logistic(x, L, k, x0=0):
    return L / (1 + np.exp(-k * (x - x0)))


def gen_logistic(p, beta_max, beta_min, center, width):
    return (beta_min + logistic(p, L=beta_max - beta_min, k=4 / width, x0=center))


def beta(p, beta_max=BETA_MAX, beta_min=BETA_MIN, center=BETA_CENTER, width=BETA_WIDTH):
    return gen_logistic(p, beta_max, beta_min, center, width)


def delta(a, delta_max=DELTA_MAX, delta_min=DELTA_MIN, center=DELTA_CENTER, width=DELTA_WIDTH):
    return gen_logistic(a, delta_max, delta_min, center, width)

# ================================================================
# 9. CALCULATE GROWTH RATES
# ================================================================

print("\n" + "=" * 70)
print("CALCULATING TRANSCRIPTIONAL GROWTH RATES")
print("=" * 70)

birth = beta(proliferation)

death = delta(apoptosis)

growth_rate = np.exp(
    birth - death
)


# Store components
adata.obs["birth_rate"] = birth

adata.obs["death_rate"] = death

adata.obs["cell_growth_rate"] = growth_rate

