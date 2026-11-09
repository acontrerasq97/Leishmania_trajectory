import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import wot
import anndata
from pathlib import Path

BASE_DIR = Path('.')

DATA_DIR = BASE_DIR / "data"
PROCESSED_DIR = DATA_DIR / "processed"
RESULTS_DIR = BASE_DIR / "results" / "wot"

TMAP_PATH = RESULTS_DIR / "tmap_leishmania"
INPUT_H5AD = PROCESSED_DIR / "ExprMatrix.embedded.h5ad"
WOT_PROBS = RESULTS_DIR / "trajectories" / "leishmania_trajectories.h5ad.txt" 
OUTPUT_FIG = BASE_DIR / "figures"/ "Para_Tiago"


adata = sc.read_h5ad(INPUT_H5AD)

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

wot_probs = pd.read_csv(WOT_PROBS, sep="\t", index_col=0)
wot_probs_norm = wot_probs.div(wot_probs.sum(axis=1), axis=0)

common_cells = adata.obs_names.intersection(wot_probs_norm.index)
adata = adata[common_cells].copy()
wot_probs_norm = wot_probs_norm.loc[common_cells]

for linaje in wot_probs_norm.columns:
    adata.obs[linaje] = wot_probs_norm[linaje]

# Ajuste global de fuentes grandes para la talk/presentación
plt.rcParams.update(
    {
        "font.size": 14,
        "axes.titlesize": 16,
        "axes.labelsize": 14,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.fontsize": 12,
        "figure.titlesize": 18,
    }
)

umap_coords = adata.obsm["X_umap"]
day_str = adata.obs["day"].astype(str)

# -------------------------------------------------------------------------
# FIGURA 1: UMAP com tipos celulares destacando as Leptomonads
# -------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 9))

temp_col = adata.obs["original_subpopulation"].astype(str).copy()

mask_not_lepto = ~temp_col.str.contains("lepto", case=False)
temp_col[mask_not_lepto] = np.nan

adata.obs["lepto_only"] = temp_col.astype("category")

lepto_colors = ["#e63946", "#457b9d", "#2a9d8f"]
palette_lepto = {
    cat: color 
    for cat, color in zip(adata.obs["lepto_only"].cat.categories, lepto_colors)
}

# 3. Graficar usando na_color para el fondo gris
sc.pl.umap(
    adata,
    color="lepto_only",
    palette=palette_lepto,
    na_color="#adb5bd",  # Color gris para los estados terminales (NaN)
    title="Leptomonads Subpopulations",
    size=35,
    legend_loc="right margin",
    na_in_legend=False,
    ax=ax,
    show=False,
)
ax.set_aspect("equal")
ax.set_xlabel("UMAP 1")
ax.set_ylabel("UMAP 2")

plt.tight_layout()
plt.savefig(OUTPUT_FIG / "1_umap_subpopulations.png", dpi=300, bbox_inches="tight")
plt.close()

# -------------------------------------------------------------------------
# FIGURA 2: Qual subpopulação de Lepto é provável ancestral (Violin / Boxplot)
# -------------------------------------------------------------------------
# Filtrar apenas as leptomonads do Dia 5
mask_lepto_d5 = (day_str == "5") & (
    adata.obs["original_subpopulation"].str.contains("lepto", case=False)
)
df_lepto_d5 = adata.obs.loc[
    mask_lepto_d5,
    ["original_subpopulation", "late_metacyclics", "haptomonads"],
].copy()
df_lepto_d5["original_subpopulation"] = df_lepto_d5[
    "original_subpopulation"
].astype(str)

# Transformar para formato longo para o Seaborn
df_melted = df_lepto_d5.melt(
    id_vars="original_subpopulation",
    value_vars=["late_metacyclics", "haptomonads"],
    var_name="Terminal_Fate",
    value_name="Probability",
)
df_melted["Terminal_Fate"] = df_melted["Terminal_Fate"].replace(
    {"late_metacyclics": "Late-Metacyclics", "haptomonads": "Haptomonads"}
)

fig, ax = plt.subplots(figsize=(10, 6))
sns.boxplot(
    data=df_melted,
    x="original_subpopulation",
    y="Probability",
    hue="Terminal_Fate",
    palette={"Late-Metacyclics": "#3a86ff", "Haptomonads": "#ff006e"},
    showmeans=True,
    meanprops={
        "marker": "o",
        "markerfacecolor": "white",
        "markeredgecolor": "black",
        "markersize": "7",
    },
    ax=ax,
)
ax.set_title("Day 5: Leptomonad Fate Commitment to Terminal States", pad=15)
ax.set_xlabel("Leptomonad Subpopulation")
ax.set_ylabel("Fate Probability")
ax.set_ylim(-0.02, 1.02)
ax.legend(title="Target Fate", frameon=True)
sns.despine()

plt.tight_layout()
plt.savefig(OUTPUT_FIG / "2_lepto_ancestor_probabilities.png", dpi=300, bbox_inches="tight")
plt.close()

# -------------------------------------------------------------------------
# FIGURA 3: Fate probabilities focadas APENAS no Dia 5 (Sem terminais)
# -------------------------------------------------------------------------
# Células do Dia 5 excluindo os próprios terminais
mask_d5_non_terminal = (
    (day_str == "5")
    & (adata.obs["original_subpopulation"] != "late_metacyclics")
    & (adata.obs["original_subpopulation"] != "haptomonads")
)

prob_meta = adata.obs["late_metacyclics"].values
prob_hapto = adata.obs["haptomonads"].values

fig, axes = plt.subplots(1, 2, figsize=(18, 8))

for idx, (vals, title, ax) in enumerate(
    [
        (prob_meta, "Day 5: Fate Probability to Late-Metacyclics", axes[0]),
        (prob_hapto, "Day 5: Fate Probability to Haptomonads", axes[1]),
    ]
):
    # Fundo de todas as células em cinza claro
    ax.scatter(
        umap_coords[:, 0],
        umap_coords[:, 1],
        c="#e9ecef",
        s=18,
        alpha=0.4,
        rasterized=True,
    )
    
    x_d5 = umap_coords[mask_d5_non_terminal, 0]
    y_d5 = umap_coords[mask_d5_non_terminal, 1]
    c_d5 = vals[mask_d5_non_terminal]

    sort_idx = np.argsort(c_d5)
    x_d5_sorted = x_d5[sort_idx]
    y_d5_sorted = y_d5[sort_idx]
    c_d5_sorted = c_d5[sort_idx]

    # Plotar as células do Dia 5 com escala completa 0 a 1
    sc_plot = ax.scatter(
        x_d5_sorted,
        y_d5_sorted,
        c=c_d5_sorted,
        cmap="coolwarm",
        s=30,
        vmin=0.0,
        vmax=1.0,
        edgecolors="none",
        alpha=0.85,
        rasterized=True,
    )
    ax.set_aspect("equal")
    ax.set_title(title, pad=12)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")

    cbar = plt.colorbar(sc_plot, ax=ax, fraction=0.035, pad=0.03)
    cbar.set_label("Probability")

plt.tight_layout()
plt.savefig(OUTPUT_FIG / "3_day5_fate_probabilities_0_to_1_coolwarm.png", dpi=300, bbox_inches="tight")
plt.close()

# -------------------------------------------------------------------------
# FIGURA 4: Log-Odds Ratio (Lineage Bias) no Dia 5 (Vermelho - Azul)
# -------------------------------------------------------------------------
# log(p_meta / p_hapto): Positivo = viés para Meta | Negativo = viés para Hapto
eps = 1e-6
log_odds = np.log2((prob_meta + eps) / (prob_hapto + eps))

# Limitar a escala entre [-3, 3] para que outliers não estourem o gradiente
limit_val = 3.0
log_odds_clipped = np.clip(log_odds, -limit_val, limit_val)

fig, ax = plt.subplots(figsize=(10.5, 9))

# Fundo cinza claro
ax.scatter(
    umap_coords[:, 0],
    umap_coords[:, 1],
    c="#e9ecef",
    s=20,
    alpha=0.4,
    rasterized=True,
)

# Células do Dia 5 com paleta divergente (Coolwarm / RdBu_r)
# Azul = Late-Metacyclics | Vermelho = Haptomonads
sc_bias = ax.scatter(
    umap_coords[mask_d5_non_terminal, 0],
    umap_coords[mask_d5_non_terminal, 1],
    c=log_odds_clipped[mask_d5_non_terminal],
    cmap="RdBu_r",
    s=30,
    vmin=-limit_val,
    vmax=limit_val,
    edgecolors="black",
    linewidths=0.25,
    rasterized=True,
)

ax.set_aspect("equal")
ax.set_title("Day 5 Lineage Bias: $\log_2(P_{meta} / P_{hapto})$", pad=15)
ax.set_xlabel("UMAP 1")
ax.set_ylabel("UMAP 2")

cbar = plt.colorbar(sc_bias, ax=ax, fraction=0.035, pad=0.03)
cbar.set_label("$\longleftarrow$ Haptomonads | Late-Metacyclics $\longrightarrow$")
cbar.set_ticks([-3, -1.5, 0, 1.5, 3])
cbar.set_ticklabels([r"$\leq -3$", "-1.5", "0", "+1.5", r"$\geq +3$"])

plt.tight_layout()
plt.savefig(OUTPUT_FIG / "4_day5_log_odds_bias_talk.png", dpi=300, bbox_inches="tight")
plt.close()


# -------------------------------------------------------------------------
# FIGURA 5: Distribuição de Ancestrais (Retrospectiva: Dia 9 -> Dias 5 e 2)
# -------------------------------------------------------------------------
targets = {"Late-Metacyclics": (day_str == "9") & (adata.obs["original_subpopulation"] == "late_metacyclics")}

anndata.read = anndata.read_h5ad

tmap_model = wot.tmap.TransportMapModel.from_directory(str(TMAP_PATH))

target_mask = targets.items()
target_cell = adata.obs_names[target_mask].tolist()   




for fate_name, target_mask in targets.items():
    target_cells = adata.obs_names[target_mask].tolist()
    cell_set_dict = {fate_name: target_cells}
    target_pop = tmap_model.population_from_cell_sets(cell_set_dict, at_time=9.0)
    
    traj_dataset = tmap_model.trajectories(target_pop)
    traj_df = traj_dataset.to_df()
    traj_prob = traj_df.sum(axis=1)

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    
    for idx, d in enumerate(["5", "2"]):
        ax = axes[idx]
        mask_day = day_str == d
        
        # 1. Capa base: TODO el UMAP del día en gris claro (esto incluye los ceros)
        ax.scatter(
            umap_coords[:, 0],
            umap_coords[:, 1],
            c="#e9ecef",
            s=20,
            alpha=0.4,
            rasterized=True,
        )
        
                # 3. FILTRAR: Solo conservar las células con probabilidad > 0
        mask_greater_than_zero = wot_probs_norm > 0
        x_filtered = umap_coords[mask_day, 0][mask_greater_than_zero]
        y_filtered = umap_coords[mask_day, 1][mask_greater_than_zero]
        p_filtered = p_day_norm[mask_greater_than_zero]
            
        # 4. ORDENAR: De menor a mayor para que los amarillos vayan al final (arriba)
        sort_idx = np.argsort(p_filtered)
        x_sorted = x_filtered[sort_idx]
        y_sorted = y_filtered[sort_idx]
        p_sorted = p_filtered[sort_idx]
        
        # 5. Graficar solo las células filtradas y ordenadas
        sc_traj = ax.scatter(
            x_sorted,
            y_sorted,
            c=p_sorted,
            cmap="viridis",     # Cambiado a viridis como en tus últimas imágenes
            s=30,
            vmin=0.0,
            vmax=1.0,           
            edgecolors="none",  
            alpha=0.9,          # Casi opaco para resaltar sobre el fondo gris
            rasterized=True,
        )
        
        ax.set_aspect("equal")
        ax.set_title(f"Day {d} Ancestors of {fate_name}")
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        
        cbar = plt.colorbar(sc_traj, ax=ax, fraction=0.035, pad=0.03)
        cbar.set_label("Relative Ancestor Probability")
        
    plt.tight_layout()
    plt.savefig(
        f"5_ancestors_{fate_name.lower().replace('-', '_')}_talk.png", 
        dpi=300, 
        bbox_inches="tight"
    )
    plt.close()

# -------------------------------------------------------------------------
# FIGURA 6: Leptomonads markers
# -------------------------------------------------------------------------

# 1. Filtrar el objeto AnnData quedándote ÚNICAMENTE con los 3 clusters de leptomonads
lepto_clusters = ["leptomonads_0", "leptomonads_1", "leptomonads_2"]
# (Asegúrate de usar la columna donde tienes estos 3 nombres, ej. 'original_subpopulation')
lepto_col = "original_subpopulation"

adata_lepto = adata[adata.obs[lepto_col].isin(lepto_clusters)].copy()
adata_lepto.obs[lepto_col] = adata_lepto.obs[lepto_col].cat.remove_unused_categories()

h5ad_file = PROCESSED_DIR / "leptomonads_subset.h5ad"
adata_lepto.write_h5ad(h5ad_file, compression="gzip")

# 2. Calcular los marcadores diferenciales (equivalente a FindAllMarkers de Seurat)
# Usando Wilcoxon rank-sum test
sc.tl.rank_genes_groups(adata_lepto, groupby=lepto_col, method="wilcoxon")

# 3. Extraer los resultados en una tabla estructurada
result = adata_lepto.uns["rank_genes_groups"]
groups = result["names"].dtype.names

markers_list = []
for group in groups:
    df_group = pd.DataFrame({
        "gene": result["names"][group],
        "logfoldchange": result["logfoldchanges"][group],
        "pvals": result["pvals"][group],
        "pvals_adj": result["pvals_adj"][group],
        "scores": result["scores"][group],
        "cluster": group
    })
    # Filtrar genes significativos (p_adj < 0.05 y sobreexpresados)
    df_sig = df_group[(df_group["pvals_adj"] < 0.05) & (df_group["logfoldchange"] > 0.25)]
    markers_list.append(df_sig)

all_markers = pd.concat(markers_list, axis=0)

# 4. Guardar en Excel y CSV
all_markers.to_csv("leptomonad_subclusters_markers.csv", index=False)
all_markers.to_excel("leptomonad_subclusters_markers.xlsx", index=False)

print("Marcadores exportados exitosamente:")
print(all_markers.groupby("cluster")["gene"].count())


import scanpy as sc

# 1. Extrair os top 5 genes de cada subcluster (mantém a ordem de importância se o df estiver ordenado)
top_markers = all_markers.groupby("cluster").head(5)["gene"].tolist()

# 2. Filtrar o AnnData para mostrar apenas as leptomonas no eixo Y (opcional, para focar a visualização)
# Se quiser mostrar todas as células, basta usar `adata` direto no dotplot
adata_lepto = adata[adata.obs["original_subpopulation"].str.contains("lepto", na=False, case=False)]

# 3. Gerar o Dotplot
sc.pl.dotplot(
    adata_lepto,
    var_names=top_markers,
    groupby="original_subpopulation",
    cmap="Reds",                     # Paleta de cores para expressão
    standard_scale="var",            # Normaliza a expressão de 0 a 1 por gene (essencial para comparar)
    title="Top Marker Genes: Leptomonad Subclusters",
    show=False
)

import matplotlib.pyplot as plt
plt.tight_layout()
plt.savefig("6_dotplot_leptomonads_markers_talk.png", dpi=300, bbox_inches="tight")
plt.show()

# 1. Isolar a lista exata dos 4 marcadores do leptomonads_2
genes_lepto2 = all_markers[all_markers["cluster"] == "leptomonads_2"]["gene"].tolist()

# 2. Plotar a expressão desses 4 genes no UMAP
sc.pl.umap(
    adata,
    color=genes_lepto2,
    cmap="magma",          # Usando a mesma paleta de probabilidades para consistência visual
    vmax="p99",            # Satura o top 1% para evitar que outliers apaguem o sinal
    size=20,
    ncols=2,               # Coloca 2 gráficos por linha (matriz 2x2)
    show=False
)

plt.suptitle("Specific Markers for leptomonads_2", fontsize=16, y=1.02)
plt.savefig("7_umap_lepto2_specific_markers_talk.png", dpi=300, bbox_inches="tight")
plt.show()