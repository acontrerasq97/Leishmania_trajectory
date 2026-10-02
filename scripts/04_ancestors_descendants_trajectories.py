
from pathlib import Path

import anndata
import pandas as pd
import wot

# ================================================================
# CONFIGURATION
# ================================================================

BASE_DIR = Path("./leishmania_trajectory")

DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"
WOT_DIR = RESULTS_DIR / "wot"

TMAP_PATH = WOT_DIR / "tmap_leishmania"

CELL_SETS_PATH = DATA_DIR / "metadata" / "cell_sets.gmt"
CELL_METADATA_PATH = DATA_DIR / "processed" / "cell_metadata.csv"
EXPRESSION_PATH = DATA_DIR / "processed" / "ExprMatrix.preprocessed.h5ad"
EXPRESSION_VAR_PATH = (DATA_DIR / "processed" / "ExprMatrix.embedded.h5ad")

adata_var = anndata.read_h5ad(EXPRESSION_VAR_PATH)

OUTPUT_DIR = WOT_DIR / "trajectories"
TRENDS_DIR = WOT_DIR / "expression_trends"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
TRENDS_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# WOT compatibility with current AnnData
# ============================================================

anndata.read = anndata.read_h5ad


# ============================================================
# Load transport maps
# ============================================================

tmap_model = wot.tmap.TransportMapModel.from_directory(str(TMAP_PATH))


# ============================================================
# Load cell metadata
# ============================================================

cell_metadata = pd.read_csv(CELL_METADATA_PATH, index_col=0)

# ================================================================
# Load expression data
# ================================================================

adata = anndata.read_h5ad(EXPRESSION_PATH)

# ============================================================
# Compatibility fix for WOT's group_cell_sets()
# ============================================================
#
# This function reproduces the logic of WOT's original
# group_cell_sets() while addressing a compatibility issue
# between the current WOT implementation and AnnData.
#
# When selecting a single cell set with:
#
#     cell_set_ds[:, i].X
#
# the resulting matrix has shape (n_cells, 1). WOT's original
# implementation uses this 2D Boolean array directly to index
# the 1D array of cell IDs, which raises an IndexError.
#
# We therefore use .flatten() to convert the Boolean mask from
# shape (n_cells, 1) to (n_cells,).
#
# The function otherwise preserves WOT's original logic:
#   1. Read the cell sets from the GMT file.
#   2. Identify the cells belonging to each cell set.
#   3. Match those cells with the verified cell metadata.
#   4. Group the cells according to their experimental day.
#   5. Store each cell set together with its corresponding time.
#
# No biological transformation or methodological modification
# is introduced. This function only fixes the indexing
# incompatibility and preserves the temporal organization
# required by WOT.
# ============================================================

def group_cell_sets_fixed(cell_set_paths, group_by_df, group_by_key="day"):

    group_to_cell_sets = {}

    if isinstance(cell_set_paths, str):
        cell_set_paths = [cell_set_paths]

    for path in cell_set_paths:

        cell_set_ds = wot.io.read_sets(path)

        for i in range(cell_set_ds.shape[1]):

            cell_set_name = cell_set_ds.var.index.values[i]

            mask = (cell_set_ds[:, i].X > 0).flatten()

            cell_ids_in_set = cell_set_ds.obs.index.values[mask]

            grouped = group_by_df[
                group_by_df.index.isin(cell_ids_in_set)
            ].groupby(group_by_key)

            for name, group in grouped:

                cell_sets = group_to_cell_sets.get(name)

                if cell_sets is None:
                    cell_sets = []
                    group_to_cell_sets[name] = cell_sets

                full_name = cell_set_name + "_" + str(name)

                cell_sets.append({
                    "set": set(group.index.values),
                    "name": full_name
                })

    return group_to_cell_sets

# ============================================================
# Group cell sets by time (day)
# ============================================================

cell_sets_by_time = group_cell_sets_fixed(str(CELL_SETS_PATH), cell_metadata, group_by_key="day")


OUTPUT_FILE = OUTPUT_DIR / "cell_sets_by_time.txt"

with open(OUTPUT_FILE, "w") as f:

    for day in sorted(cell_sets_by_time.keys()):

        f.write(f"\nDay {day}\n")

        for s in cell_sets_by_time[day]:

            f.writelines(f'{s["name"]}\t{len(s["set"])}\n')


# ============================================================
# Calculate trajectories and expression trends
# ============================================================

for day, cell_sets in sorted(cell_sets_by_time.items()):

    # Convert cell sets to the dictionary format expected by WOT
    cell_sets = {
        s["name"]: list(s["set"])
        for s in cell_sets_by_time[day]
    }

    # Create population vectors
    populations = tmap_model.population_from_cell_sets(cell_sets, at_time=day)

    # Calculate trajectories
    trajectory_lmj = tmap_model.trajectories(populations)

    # Save trajectories
    #output_path = OUTPUT_DIR / f"leishmania_trajectories_day_{day}.h5ad"

    #trajectory_lmj.write_h5ad(output_path)

    # Calculate expression trends along trajectories
    trajectory_trends = (wot.tmap.trajectory_trends_from_trajectory(trajectory_lmj, adata))

    # Save expression trends
    #for i, trend in enumerate(trajectory_trends):

    #    trends_path = (
    #        TRENDS_DIR
    #        / f"expression_trends_day_{day}_{i}.h5ad"
    #    )

    #    trend.write_h5ad(trends_path)

# ============================================================
# Heatmap of expression trends
# ============================================================

import matplotlib.pyplot as plt
import seaborn as sns

trend = trajectory_trends[0]

# Convert expression matrix to DataFrame
trend_df = pd.DataFrame(
    trend.X,
    index=trend.obs_names,
    columns=trend.var_names
)

# Transpose:
# rows    = genes
# columns = days
heatmap_df = trend_df.T

# Select genes with the largest variation across time
gene_variation = heatmap_df.max(axis=1) - heatmap_df.min(axis=1)

top_genes = gene_variation.sort_values(
    ascending=False
).head(100).index

heatmap_df = heatmap_df.loc[top_genes]

# Plot
plt.figure(figsize=(8, 14))

sns.heatmap(
    heatmap_df,
    cmap="viridis",
    xticklabels=True,
    yticklabels=True
)

plt.xlabel("Day")
plt.ylabel("WOT-weighted log-expression")
plt.title("Expression trends along WOT trajectories")

plt.tight_layout()

FIGURE_PATH = TRENDS_DIR / "expression_trends_heatmap_top50.png"

plt.savefig(
    FIGURE_PATH,
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print(f"Heatmap saved to: {FIGURE_PATH}")


# ============================================================
# Shared Ancestry / Trajectory divergence
# ============================================================

adata_var = anndata.read_h5ad(EXPRESSION_VAR_PATH)
