
from pathlib import Path

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import anndata
import wot
import seaborn as sns


#===================================================
#CONFIGURATION
#===================================================

BASE_DIR = Path("./leishmania_trajectory")

WOT_DIR = BASE_DIR / "results"/ "wot"
METADATA_DIR = BASE_DIR / "data" / "metadata"
PROCESSED_DIR = BASE_DIR / "data"/ "processed"

TMAP_PATH = WOT_DIR / "tmap_leishmania"
CELL_SETS_PATH = METADATA_DIR / "cell_sets.gmt"

CELL_METADATA_PATH = PROCESSED_DIR / "cell_metadata.csv"
OUTPUT_DIR = WOT_DIR / "fate_matrices"
FATE_H5AD_PATH = OUTPUT_DIR / "fate_matrices_model.h5ad"
FATE_CSV_PATH = OUTPUT_DIR / "fate_matrices_table.csv"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

#====================================================
# WOT COMPATIBILITY
#====================================================
anndata.read = anndata.read_h5ad


# Load transport map model and cell sets
tmap_model = wot.tmap.TransportMapModel.from_directory(TMAP_PATH)
cell_sets = wot.io.read_sets(CELL_SETS_PATH, as_dict=True)

# create indicator vectors for each cell set
target_destinations = tmap_model.population_from_cell_sets(cell_sets, at_time=9)

# ============================================================
# CALCULATE FATE MATRICES
# ============================================================

fate_lmj = tmap_model.fates(target_destinations)

fate_lmj.write_h5ad(FATE_H5AD_PATH)

fate_df = fate_lmj.to_df() if hasattr(fate_lmj, "to_df") else pd.DataFrame(fate_lmj.X, index=fate_lmj.obs_names, columns=fate_lmj.var_names)

fate_df.to_csv(FATE_CSV_PATH)
