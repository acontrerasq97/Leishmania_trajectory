import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

from pathlib import Path


BASE_DIR = Path('./')

DATA_DIR = BASE_DIR / "data"
PROCESSED_DIR = DATA_DIR / "processed"
RESULTS_DIR = BASE_DIR / "results" / "wot"

INPUT_H5AD = PROCESSED_DIR / "ExprMatrix.embedded.h5ad"
WOT_PROBS = RESULTS_DIR / "trajectories" / "leishmania_trajectories.h5ad.txt" 


adata = sc.read_h5ad(INPUT_H5AD)
wot_probs = pd.read_csv(WOT_PROBS, sep="\t", index_col=0)

common_cells = adata.obs_names.intersection(wot_probs.index)
adata = adata[common_cells].copy()
wot_probs = wot_probs.loc[common_cells]

for linaje in wot_probs.columns:
    adata.obs[linaje] = wot_probs[linaje]

tiempo_wot = adata.obs["day"].astype(float)

adata.obs["dominant_dest"] = wot_probs.idxmax(axis=1)
sc.pp.neighbors(adata, n_neighbors=30, use_rep="X_pca")
sc.tl.leiden(adata, resolution=1.0, key_added="subestados")
sc.tl.paga(adata, groups="subestados")

adata.obs["day_numeric"] = adata.obs["day"].astype(float)
sc.pl.paga(
    adata,
    color=["day", "destino_dominante"],
    title=["Day", "Sesgo de Linaje WOT"],
    threshold=0.15,  # Umbral de conectividad entre ramas
    layout="fr",      # Fruchterman-Reingold (layout dirigido por fuerzas similar al de URD)
    show=False
)

sc.tl.draw_graph(adata, init_pos="paga")
