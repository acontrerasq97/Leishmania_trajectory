from pathlib import Path

import pandas as pd
import scanpy as sc
import numpy as np


BASE_DIR = Path("./leishmania_trajectory")

DATA_DIR = BASE_DIR / 'data'
H5AD_DIR = DATA_DIR / 'processed'/ 'ExprMatrix.embedded.h5ad'
META_DIR = DATA_DIR / 'metadata' / 'cell_sets.gmt'

adata = sc.read_h5ad(H5AD_DIR)

adata.obs['cell_type'] = np.nan

with open(META_DIR, 'r') as f:
    for line in f:
        
        parts = line.strip().split('\t')
        if len(parts) < 3:
            continue
            
        group_name = parts[0]   
        cell_ids = parts[2:]   
        
        adata.obs.loc[adata.obs.index.isin(cell_ids), 'cell_type'] = group_name

# 4. Convertir a categoría para que Scanpy lo procese eficientemente
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

#Wilcoxon test
sc.tl.rank_genes_groups(
    adata, 
    groupby='cell_type', 
    reference='leptomonads_1', 
    method='wilcoxon',
    key_added='rank_leptomonads_2'
)

result_df = sc.get.rank_genes_groups_df(adata, group="leptomonads_2", key="rank_leptomonads_2")

up_regulated_genes = result_df[
    (result_df['logfoldchanges'] > 0.25) & 
    (result_df['pvals_adj'] < 0.01)
].copy()

top_genes =  up_regulated_genes.sort_values(by='logfoldchanges', ascending=False )

print(result_df.head(10)[['names', 'logfoldchanges', 'pvals_adj']])
