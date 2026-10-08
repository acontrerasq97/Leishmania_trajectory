
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
TMAP_PATH = WOT_DIR / "tmap_leishmania"

RBPS_OUTPUT = WOT_DIR / 'predictive_rbps'
GENES_OUTPUT = WOT_DIR / 'predictive_genes'
RBPS_OUTPUT.mkdir(parents=True, exist_ok=True)
GENES_OUTPUT.mkdir(parents=True, exist_ok=True)

anndata.read = anndata.read_h5ad

# Load expression dataset and subset to transcription factors
adata = wot.io.read_dataset(FULL_DS_PATH, var_filter=GENES_PATH)

# Load transport map model and cell sets
tmap_model = wot.tmap.TransportMapModel.from_directory(TMAP_PATH)
cell_sets = wot.io.read_sets(CELL_SETS_PATH, as_dict=True)

# create indicator vector for haptomonads cell set at day 9
target_cell_set = tmap_model.population_from_cell_sets({'late_metacyclics':cell_sets['late_metacyclics']}, at_time=9)

# Compute fate matrix for haptomonads 
fate_ds = tmap_model.fates(target_cell_set)

# Find differentially expressed genes at day 5
results = wot.tmap.diff_exp(adata[adata.obs['day'].isin([5])], fate_ds, compare='all')
results[(results['t_fdr']<0.01)&(results['name1']=='late_metacyclics')].sort_values('fraction_expressed_ratio', ascending=False)


df_results = results.copy()
df_map = pd.read_csv(DATA_DIR / 'metadata' / 'genes_symbol.csv')
df_map = df_map[['Gene ID', 'Gene Name or Symbol', 'Product Description']].copy()
df_map.columns = ['gene_id', 'gene_symbol', 'description']

df_final = pd.merge(df_results, df_map, left_index=True, right_on='column00', how='left')

#results.to_csv(RBPS_OUTPUT / 'RBPS_haptomonads_9_5.csv')
df_final.to_csv(GENES_OUTPUT / 'genes_late_metacyclics_9_5.csv')


import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Set visual theme for scientific publication
sns.set_theme(style="ticks", context="talk")
plt.figure(figsize=(11, 8))

# Calculate -log10 of the FDR/adjusted p-value
# A small delta (1e-30) avoids log(0) errors if t_fdr is exactly 0
df_final['-log10_fdr'] = -np.log10(df_final['t_fdr'] + 1e-30)

# Define standard statistical thresholds for coloring
threshold_fdr = -np.log10(0.05) # Equivalent to FDR < 0.05
threshold_fc_up = 1.05
threshold_fc_down = 0.95

# Color code points based on biological regulation status
colors = []
for idx, row in df_final.iterrows():
    if row['t_fdr'] < 0.05 and row['fold_change'] >= threshold_fc_up:
        colors.append('#E63946') # Upregulated (Red)
    elif row['t_fdr'] < 0.05 and row['fold_change'] <= threshold_fc_down:
        colors.append('#457B9D') # Downregulated (Blue)
    else:
        colors.append('#A8DADC') # Not Significant (Gray)

# Draw the Volcano scatter plot
plt.scatter(df_final['fold_change'], df_final['-log10_fdr'], 
            c=colors, alpha=0.8, edgecolors='none', s=60)

# Add dashed threshold guide lines
plt.axhline(y=threshold_fdr, color='gray', linestyle='--', linewidth=1.2)
plt.axvline(x=threshold_fc_up, color='gray', linestyle='--', linewidth=1.2)
plt.axvline(x=threshold_fc_down, color='gray', linestyle='--', linewidth=1.2)

# --- GENE LABELING LOGIC FOR BOTH UP AND DOWNREGULATED ---
# Filter for statistically significant entries (FDR < 0.05)
significant_genes = df_final[df_final['t_fdr'] < 0.05]

# Select top 3 highest Fold Change (Upregulated)
top_up = significant_genes.sort_values(by='fold_change', ascending=False).head(3)

# Select top 3 lowest Fold Change (Downregulated)
top_down = significant_genes.sort_values(by='fold_change', ascending=True).head(3)

# Combine both sets to apply text labels on the plot
genes_to_label = pd.concat([top_up, top_down])

for idx, row in genes_to_label.iterrows():
    # Fallback to genomic LmjF ID if gene symbol is missing or N/A
    label = row['gene_symbol'] if str(row['gene_symbol']) != 'N/A' and pd.notna(row['gene_symbol']) else row['gene_id']
    
    # Adjust horizontal text alignment based on which side of the plot it falls on
    offset_x = 0.002 if row['fold_change'] >= 1 else -0.015
    
    plt.text(row['fold_change'] + offset_x, row['-log10_fdr'] + 0.2, label, 
             fontsize=11, weight='bold', alpha=0.85, 
             ha='left' if row['fold_change'] >= 1 else 'right')

# Axis and Title Customization
plt.title('Differentially Regulated RBPs: Haptomonads vs Other', pad=15, weight='bold')
plt.xlabel('Fold Change (Relative Expression)')
plt.ylabel(r'Significance: $-log_{10}(FDR)$')

# Retrieve current axis limits to place descriptive annotations in the top corners
x_min, x_max = plt.xlim()
y_max = plt.ylim()[1]

# Add corner text indicators for clarity
plt.text(x_min + 0.005, y_max - 1.5, "⬇️ Downregulated\n(Turned off in trajectory)", color='#457B9D', weight='bold', fontsize=12)
plt.text(x_max - 0.005, y_max - 1.5, "⬆️ Upregulated\n(Turned on in trajectory)", color='#E63946', weight='bold', fontsize=12, ha='right')

sns.despine()
plt.tight_layout()

# Save high-resolution figure for publication
plt.savefig('volcano_plot_rbps_english.png', dpi=300)
plt.show()
