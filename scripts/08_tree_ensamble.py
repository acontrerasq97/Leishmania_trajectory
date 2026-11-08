import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from pathlib import Path


BASE_DIR = Path('./leishmania_trajectory')

DATA_DIR = BASE_DIR / "data"
PROCESSED_DIR = DATA_DIR / "processed"
RESULTS_DIR = BASE_DIR / "results" / "wot"

INPUT_H5AD = PROCESSED_DIR / "ExprMatrix.embedded.h5ad"
WOT_PROBS = RESULTS_DIR / "trajectories" / "leishmania_trajectories.h5ad.txt" 
OUTPUT_FIG = BASE_DIR / "figures"/ "wot"


adata = sc.read_h5ad(INPUT_H5AD)
wot_probs = pd.read_csv(WOT_PROBS, sep="\t", index_col=0)
wot_probs_norm = wot_probs.div(wot_probs.sum(axis=1), axis=0)

common_cells = adata.obs_names.intersection(wot_probs_norm.index)
adata = adata[common_cells].copy()
wot_probs_norm = wot_probs_norm.loc[common_cells]

for linaje in wot_probs_norm.columns:
    adata.obs[linaje] = wot_probs_norm[linaje]

adata.obs["day"] = adata.obs["day"].astype("category")

cmap_ember = LinearSegmentedColormap.from_list(
    "ember_clean", ["#030101", "#f6d365", "#fda085", "#d63031"]
)

# Paleta personalizada: Gris neutro -> Turquesa -> Azul cobalto -> Violeta
cmap_neon_blue = LinearSegmentedColormap.from_list(
    "neon_blue", ["#020304", "#00cec9", "#0984e3", "#8e82e8"]
)

prob_threshold = 0.5

probs_haptomonads = adata.obs["haptomonads"].copy()
probs_haptomonads[probs_haptomonads <= prob_threshold] = np.nan
adata.obs["haptomonads_threshold"] = probs_haptomonads

probs_metacyclics = adata.obs["late_metacyclics"].copy()
probs_metacyclics[probs_metacyclics <= prob_threshold] = np.nan
adata.obs["late_metacyclics_threshold"] = probs_metacyclics

fig, axes = plt.subplots(1, 3, figsize=(14, 6))

sc.pl.umap(
    adata,
    color="day",
    title="Sampling Day",
    size=10,
    ax=axes[0],
    show=False
)
axes[0].set_xlabel("UMAP 1")
axes[0].set_ylabel("UMAP 2")

sc.pl.umap(
    adata,
    color="late_metacyclics_threshold",
    title="Fate Probability towards Late-Metacyclics",
    cmap="magma",        
    na_color="#ececec",
    size=15,
    sort_order=True,        
    vmax="p99",             
    ax=axes[1],
    show=False
)
axes[1].set_xlabel("UMAP 1")
axes[1].set_ylabel("UMAP 2")

sc.pl.umap(
    adata,
    color="haptomonads_threshold",
    title="Fate Probability towards Haptomonads",
    cmap="magma",       
    na_color="#ececec",   
    size=15,
    sort_order=True,        
    vmax="p99",             
    ax=axes[2],
    show=False
)
axes[2].set_xlabel("UMAP 1")
axes[2].set_ylabel("UMAP 2")

plt.tight_layout()

plt.savefig(
    OUTPUT_FIG / "leishmania_trajectories_figure.pdf",
    format="pdf",
    bbox_inches="tight",  
    transparent=False,
)
plt.savefig(
    OUTPUT_FIG / "leishmania_trajectories_figure_4K.png",
    dpi=600,  
    bbox_inches="tight",
    facecolor="white",  
    transparent=False,
)

plt.show()



# -------------------------------------------------------------
# 1. Parámetros y criterios comunes
# -------------------------------------------------------------
# Diccionario para almacenar el mapeo: {código_de_barra: nombre_subpoblacion}
cell_to_subpop = {}

# Lectura del archivo GMT/GTM
with open(DATA_DIR / "metadata" / "cell_sets.gmt", "r") as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) >= 3:
            subpop_name = parts[0]
            # En formato GMT estándar: col 0 = grupo, col 1 = descripción, col 2 en adelante = barcodes
            barcodes = parts[2:]
            for bc in barcodes:
                cell_to_subpop[bc] = subpop_name

# Mapear directamente a adata.obs preservando el índice original de células
adata.obs["original_subpopulation"] = adata.obs_names.map(cell_to_subpop).astype("category")

# Verificar que todas las células se hayan asignado
print("Distribución de células por subpoblación original:")
print(adata.obs["original_subpopulation"].value_counts(dropna=False))

# Filtro temporal de precursores (Días 2 y 5) y umbral estricto
is_progenitor = adata.obs["day"].isin([2, 5])
threshold = 0.50

umap_x = adata.obsm["X_umap"][:, 0]
umap_y = adata.obsm["X_umap"][:, 1]

# -------------------------------------------------------------
# A. Tabla de Precursores hacia Late-Metacyclics
# -------------------------------------------------------------
mask_meta = is_progenitor & (adata.obs["late_metacyclics"] > threshold)

df_precursors_meta = pd.DataFrame(
    {
        "cell_barcode": adata.obs_names[mask_meta],
        "day": adata.obs.loc[mask_meta, "day"].values,
        "original_subpopulation": adata.obs.loc[mask_meta, "original_subpopulation"].values,
        "umap_1": umap_x[mask_meta],
        "umap_2": umap_y[mask_meta],
        "prob_late_metacyclics": adata.obs.loc[mask_meta, "late_metacyclics"].values,
        "prob_haptomonads": adata.obs.loc[mask_meta, "haptomonads"].values,
    },
    index=adata.obs_names[mask_meta],
)

# -------------------------------------------------------------
# B. Tabla de Precursores hacia Haptomonads
# -------------------------------------------------------------
mask_hapto = is_progenitor & (adata.obs["haptomonads"] > threshold)

df_precursors_hapto = pd.DataFrame(
    {
        "cell_barcode": adata.obs_names[mask_hapto],
        "day": adata.obs.loc[mask_hapto, "day"].values,
        "original_subpopulation": adata.obs.loc[mask_hapto, "original_subpopulation"].values,
        "umap_1": umap_x[mask_hapto],
        "umap_2": umap_y[mask_hapto],
        "prob_haptomonads": adata.obs.loc[mask_hapto, "haptomonads"].values,
        "prob_late_metacyclics": adata.obs.loc[mask_hapto, "late_metacyclics"].values,
    },
    index=adata.obs_names[mask_hapto],
)

# -------------------------------------------------------------
# C. Guardar archivos en formato CSV
# -------------------------------------------------------------
df_precursors_meta.to_csv(RESULTS_DIR / "precursors_late_metacyclics.csv", index=False)
df_precursors_hapto.to_csv(RESULTS_DIR / "precursors_haptomonads.csv", index=False)

print(f"\nTotal precursores Late-Metacyclics guardados: {len(df_precursors_meta)}")
print(f"Total precursores Haptomonads guardados: {len(df_precursors_hapto)}")

print("\n--- Desglose de precursores Late-Metacyclics por subpoblación original ---")
print(df_precursors_meta["original_subpopulation"].value_counts())

print("\n--- Desglose de precursores Haptomonads por subpoblación original ---")
print(df_precursors_hapto["original_subpopulation"].value_counts())
