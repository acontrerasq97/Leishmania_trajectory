

from pathlib import Path

import anndata
import wot

# ================================================================
# CONFIGURATION
# ================================================================

BASE_DIR = Path("./leishmania_trajectory")

DATA_DIR = BASE_DIR / "data" / "metadata"
RESULTS_DIR = BASE_DIR / "results"
WOT_DIR = RESULTS_DIR / "wot"

TMAP_PATH = WOT_DIR / "tmap_leishmania"

CELL_SETS_PATH = DATA_DIR / "cell_sets.gmt"

OUTPUT_DIR = WOT_DIR / "trajectories"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# Load transport maps
# ============================================================

tmap_model = wot.tmap.TransportMapModel.from_directory(TMAP_PATH)

# ============================================================
# Load cell sets
# ============================================================

cell_sets = wot.io.read_sets(CELL_SETS_PATH, as_dict=True)


# ============================================================
# Create population vectors
# ============================================================

populations = tmap_model.population_from_cell_sets(cell_sets, at_time=9)

# ============================================================
# Calculate trajectories
# ============================================================

# Compatibility with WOT
anndata.read = anndata.read_h5ad
trajectory_lm = tmap_model.trajectories(populations)

# ============================================================
# Save trajectories
# ============================================================
output_path = OUTPUT_DIR / "leishmania_trajectories.h5ad"

wot.io.write_dataset(trajectory_lm, str(output_path))
