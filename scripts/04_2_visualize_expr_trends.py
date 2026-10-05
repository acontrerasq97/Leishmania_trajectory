from pathlib import Path

from IPython.display import display
import anndata
import matplotlib.pyplot as plt
import ipywidgets as widgets


# ================================================================
# CONFIGURATION
# ================================================================

TRENDS_DIR = Path(
    "leishmania_trajectory/results/wot/expression_trends"
)


# ================================================================
# LOAD EXPRESSION TREND DATASETS
# ================================================================

trend_paths = sorted(TRENDS_DIR.glob("expression_trends_*.h5ad"))

trajectory_trend_datasets = [
    anndata.read_h5ad(path)
    for path in trend_paths
]

trajectory_names = [
    path.stem.replace("expression_trends_", "")
    for path in trend_paths
]


print(f"Number of trajectories: {len(trajectory_trend_datasets)}")
print(trajectory_names)


# ================================================================
# TRAJECTORY SELECTION
# ================================================================

trajectory_dropdown = widgets.SelectMultiple(
    options=trajectory_names,
    value=(trajectory_names[0],),
    description="Trajectory:"
)


# ================================================================
# GENE INPUT
# ================================================================

gene_input = widgets.Text(
    placeholder='',
    description="Genes:",
    value="Nanog",
    continuous_update=False
)


# ================================================================
# PLOT EXPRESSION TRENDS
# ================================================================

def update_trends_vis(selected_trajectories, gene_names):

    gene_names = (gene_names.replace(" ", "").split(","))

    gene_names_lower = [gene.lower() for gene in gene_names]

    figure = plt.figure(figsize=(10, 7))

    for selected_trajectory in selected_trajectories:
        trajectory_index = trajectory_names.index(selected_trajectory)
        mean = trajectory_trend_datasets[trajectory_index]

        # Select requested genes
        gene_mask = (mean.var_names.str.lower().isin(gene_names_lower))

        mean = mean[:, gene_mask]

        if mean.shape[1] == 0:
            continue

        # Days
        timepoints = mean.obs_names.astype(float)

        # Plot each selected gene
        for i in range(mean.shape[1]):

            mean_i = mean[:, i].X.flatten()

            plt.plot(timepoints, 
                    mean_i, 
                    marker="o", 
                    label=(mean.var_names[i] + " — " + selected_trajectory))

    plt.xlabel("Day")
    plt.ylabel("WOT-weighted log-expression")
    plt.xticks([0, 2, 5, 9])
    plt.legend()
    plt.tight_layout()

    display(fig)
    return figure

# ================================================================
# INTERACTIVE WIDGET
# ================================================================

widgets.interact(
    update_trends_vis,
    selected_trajectories=trajectory_dropdown,
    gene_names=gene_input
)
