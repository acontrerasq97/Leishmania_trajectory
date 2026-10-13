
from pathlib import Path

import anndata
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import wot

BASE_DIR = Path ('leishmania_trajectory/')
WOT_DIR = BASE_DIR / "results"/ "wot"

DATA_DIR = BASE_DIR / 'data'

FULL_DS_PATH = DATA_DIR / 'principal' / 'ExprMatrix.h5ad'
CELL_DAYS_PATH = DATA_DIR / 'metadata'/ 'cell_days.txt'
CELL_SETS_PATH = DATA_DIR / 'metadata'/ 'cell_sets.gmt'
GENES_PATH = DATA_DIR / 'metadata'/ 'genes.txt'
RBPS_PATH = DATA_DIR / 'metadata' / 'RBPs.txt'
TMAP_PATH = WOT_DIR / "tmap_leishmania"

RBPS_OUTPUT = WOT_DIR / 'predictive_rbps'
GENES_OUTPUT = WOT_DIR / 'predictive_genes'
RBPS_OUTPUT.mkdir(parents=True, exist_ok=True)
GENES_OUTPUT.mkdir(parents=True, exist_ok=True)

#====================================================
# WOT COMPATIBILITY
#====================================================
anndata.read = anndata.read_h5ad

# Load expression dataset and subset to RBPs
adata = wot.io.read_dataset(FULL_DS_PATH, var_filter=RBPS_PATH)
# Load expression dataset and subset to genes
#adata = wot.io.read_dataset(FULL_DS_PATH, var_filter=GENES_PATH)

# Load transport map model and cell sets
tmap_model = wot.tmap.TransportMapModel.from_directory(TMAP_PATH)
cell_sets = wot.io.read_sets(CELL_SETS_PATH, as_dict=True)

# create indicator vector for haptomonads cell set at day 9
target_cell_set = tmap_model.population_from_cell_sets({'haptomonads':cell_sets['haptomonads']}, at_time=9)

# Compute fate matrix for haptomonads 
fate_ds = tmap_model.fates(target_cell_set)

# Find differentially expressed genes at day 5
results = wot.tmap.diff_exp(adata[adata.obs['day'].isin([5])], fate_ds, compare='all')
results[(results['t_fdr']<0.01)&(results['name1']=='haptomonads')].sort_values('fraction_expressed_ratio', ascending=False)
results_top10 =  results[(results['t_fdr']<0.01)&(results['name1']=='haptomonads')].sort_values('fraction_expressed_ratio', ascending=False)

df_results_top10 = results_top10.copy()
df_map = pd.read_csv(DATA_DIR / 'metadata' / 'genes_symbol.csv')
df_map = df_map[['Gene ID', 'Gene Name or Symbol', 'Product Description']].copy()
df_map.columns = ['gene_id', 'gene_symbol', 'description']

df_final = pd.merge(df_results_top10, df_map, left_index=True, right_on='gene_id', how='left')

#df_final.to_csv(RBPS_OUTPUT / 'RBPS_haptomonads_9_5.csv')
#df_final.to_csv(GENES_OUTPUT / 'genes_late_metacyclics_9_5.csv')

# Set visual theme for scientific publication
sns.set_theme(style="ticks", context="talk")
plt.figure(figsize=(11, 8))

# Calculate -log10 of the FDR/adjusted p-value
# A small delta (1e-30) avoids log(0) errors if t_fdr is exactly 0
df_final['-log10_fdr'] = -np.log10(df_final['t_fdr'] + 1e-30)

# Define standard statistical thresholds for coloring
threshold_fdr = -np.log10(0.05) # Equivalent to FDR < 0.05
threshold_fc_up = 1.05
threshold_fc_down = 1

# Color code points based on biological regulation status
colors = []
for idx, row in df_final.iterrows():
    if row['t_fdr'] < 0.05 and row['fold_change'] >= threshold_fc_up:
        colors.append('#E63946') # Upregulated (Red)
    elif row['t_fdr'] < 0.05 and row['fold_change'] < threshold_fc_down:
        colors.append('#457B9D') # Downregulated (Blue)
    else:
        colors.append('#A8DADC') # Not Significant (Gray)

# Draw the Volcano scatter plot
plt.scatter(df_final['fold_change'], df_final['-log10_fdr'], 
            c=colors, alpha=0.8, edgecolors='none', s=15)

# Add dashed threshold guide lines
plt.axhline(y=threshold_fdr, color='gray', linestyle='--', linewidth=1.2)
plt.axvline(x=threshold_fc_up, color='gray', linestyle='--', linewidth=1.2)
plt.axvline(x=threshold_fc_down, color='gray', linestyle='--', linewidth=1.2)

# --- GENE LABELING LOGIC FOR BOTH UP AND DOWNREGULATED ---
# Filter for statistically significant entries (FDR < 0.05)
significant_genes = df_final[df_final['t_fdr'] < 0.05]

# Select top 15 highest Fold Change (Upregulated)
top_up = significant_genes.sort_values(by='fold_change', ascending=False).head(10)

# Select top 15 lowest Fold Change (Downregulated)
top_down = significant_genes.sort_values(by='fold_change', ascending=True).head(10)

# Combine both sets to apply text labels on the plot
genes_to_label = pd.concat([top_up, top_down])

from adjustText import adjust_text

texts = []

for idx, row in genes_to_label.iterrows():
    # Fallback to genomic LmjF ID if gene symbol is missing or N/A
    label = row['gene_symbol'] if str(row['gene_symbol']) != 'N/A' and pd.notna(row['gene_symbol']) else row['gene_id']
      
    t = plt.text(row['fold_change'], row['-log10_fdr'] + 0.2, label, 
             fontsize=8, weight='bold', alpha=0.85)
    texts.append(t)

adjust_text(texts, 
            arrowprops=dict(arrowstyle='-', color='gray', lw=0.5),
            expand=(1.2, 1.5))
# Axis and Title Customization
plt.title('Differentially Regulated Genes: Late-Metacyclics trajectory', pad=15, weight='bold')
plt.xlabel('Fold Change (Relative Expression)')
plt.ylabel(r'Significance: $-log_{10}(FDR)$')

sns.despine()
plt.tight_layout()

# Save high-resolution figure for publication
#plt.savefig(BASE_DIR / 'figures' / 'wot' / 'VP_rbps_haptomonads.png', dpi=300)
# Save high-resolution figure for publication
plt.savefig(BASE_DIR / 'figures' / 'wot' / 'VP_genes_late_metacyclics.png', dpi=300)


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import MinMaxScaler

# 1. Definir tus 10 RBPs más significativas (usando los IDs reales de tu tabla de WOT)
top_rbps = [
    "LmjF.33.1460", "LmjF.04.1170", "LmjF.18.0590", "LmjF.24.1570", 
    "LmjF.35.4130", "LmjF.07.0340", "LmjF.35.2200", "LmjF.36.1635", 
    "LmjF.32.0750", "LmjF.36.0050"
]

# 2. Extraer la expresión de tus RBPs y asociar los metadatos de las células
# Extraemos la matriz de expresión desde tu objeto AnnData de Scanpy
df_expression = adata[:, top_rbps].to_df()
df_expression['day'] = adata.obs['day'].values

# Extraer el peso de probabilidad hacia haptomónadas calculado por WOT
# Opción A (La más recomendada si es un array de NumPy):
df_expression['trajectory_weight'] = fate_ds[:, 'haptomonads'].X.flatten()

# Opción B (Si el objeto es una matriz dispersa o requiere conversión densa previa):
# df_expression['trajectory_weight'] = np.asarray(fate_ds[:, 'haptomonads'].X).squeeze()

# 3. Calcular la expresión promedio ponderada por día
# Multiplicamos la expresión de cada gen en cada célula por su peso en la trayectoria
time_trends = {}
unique_days = sorted(df_expression['day'].unique())

for gene in top_rbps:
    weighted_means = []
    for day in unique_days:
        df_day = df_expression[df_expression['day'] == day]
        # Si la suma de pesos del día es cero, evitamos división por cero
        total_weight = df_day['trajectory_weight'].sum()
        if total_weight > 0:
            weighted_mean = (df_day[gene] * df_day['trajectory_weight']).sum() / total_weight
        else:
            weighted_mean = df_day[gene].mean()
        weighted_means.append(weighted_mean)
    time_trends[gene] = weighted_means

# Crear el DataFrame de tendencias temporales (Filas = Genes, Columnas = Días)
df_trends = pd.DataFrame(time_trends, index=[f"Day {int(d)}" for d in unique_days]).T

# 4. NORMALIZACIÓN (Escala Min-Max por fila)
# Para poder comparar genes muy abundantes (como LmjF.32.0750 con media 2.5) con genes de baja expresión (como LmjF.04.1170 con media 0.16)
# llevamos todos los valores a una escala de 0 (expresión mínima) a 1 (expresión máxima del gen en el tiempo).
scaler = MinMaxScaler()
df_normalized = pd.DataFrame(
    scaler.fit_transform(df_trends.T).T, 
    index=df_trends.index, 
    columns=df_trends.columns
)

# 5. ORDENAR LOS GENES EN CASCADA TEMPORAL (El truco del pico de expresión)
# Buscamos en qué columna (día) cada gen alcanza su valor máximo (1.0)
peak_days = df_normalized.idxmax(axis=1)

# Mapeamos los nombres de los días a un valor numérico para poder ordenar jerárquicamente
day_to_num = {col: i for i, col in enumerate(df_normalized.columns)}
df_normalized['peak_order'] = peak_days.map(day_to_num)

# Ordenamos el DataFrame de menor a mayor día de pico de expresión y eliminamos la columna de orden
df_heatmap = df_normalized.sort_values(by='peak_order', ascending=True).drop(columns=['peak_order'])

# --- REEMPLAZO OPCIONAL POR NOMBRES COMUNES SI HICISTE EL MERGE ---
# Si quieres que en el gráfico se lean los nombres comunes (como PABP1) en lugar de LmjF:
# mapeo_nombres = dict(zip(df_annotation['gene_id'], df_annotation['gene_symbol']))
# df_heatmap.index = [mapeo_nombres.get(idx, idx) for idx in df_heatmap.index]

# 6. DIBUJAR EL HEATMAP
plt.figure(figsize=(10, 8))
sns.set_theme(style="white", context="talk")

# Usamos la paleta 'YlOrRd' (Amarillo-Naranja-Rojo) o 'viridis' para un look clásico de transcriptómica
ax = sns.heatmap(
    df_heatmap, 
    cmap="YlOrRd", 
    linewidths=0.5, 
    cbar_kws={'label': 'Relative Expression (Min-Max scaled)'},
    annot=False
)

plt.title('Trajectory Cascade Heatmap: RBPs Transition to Haptomonads', pad=20, weight='bold')
plt.xlabel('Developmental Time', labelpad=15)
plt.ylabel('RNA-Binding Proteins (Ordered by Peak)', labelpad=15)

plt.tight_layout()
plt.savefig('trajectory_heatmap_cascade.png', dpi=300)
plt.show()
