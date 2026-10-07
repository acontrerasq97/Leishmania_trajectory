
from pathlib import Path

import anndata
import pandas as pd
from pybiomart import Server
import wot

BASE_DIR = Path ('leishmania_trajectory/')
WOT_DIR = BASE_DIR / "results"/ "wot"

DATA_DIR = BASE_DIR / 'data'

FULL_DS_PATH = DATA_DIR / 'principal' / 'ExprMatrix.h5ad'
CELL_DAYS_PATH = DATA_DIR / 'metadata'/ 'cell_days.txt'
CELL_SETS_PATH = DATA_DIR / 'metadata'/ 'cell_sets.gmt'
RBPS_PATH = DATA_DIR / 'metadata'/ 'RBPs.txt'
TMAP_PATH = WOT_DIR / "tmap_leishmania"

RBPS_OUTPUT = WOT_DIR / 'predictive_rbps'
RBPS_OUTPUT.mkdir(parents=True, exist_ok=True)

anndata.read = anndata.read_h5ad

# Load expression dataset and subset to transcription factors
adata = wot.io.read_dataset(FULL_DS_PATH, var_filter=RBPS_PATH)

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

results.to_csv(RBPS_OUTPUT / 'RBPS_haptomonads_9_5.csv')


server = Server(host='protists.ensembl.org')
dataset = server.marts['protists_mart'].datasets['lmajor_eg_gene']

print(server.marts['protists_mart'].datasets.keys())
queries = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name', 'description'])

