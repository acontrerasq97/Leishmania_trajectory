
import ipywidgets as widgets
import matplotlib.pyplot as plt
import numpy as np
import math

# 1. Configuración de los menús interactivos vinculados a tu objeto 'fate_lmj'
fate_dropdown1 = widgets.Dropdown(
    options=fate_lmj.var.index,
    description='Fate 1:'
)
fate_dropdown2 = widgets.Dropdown(
    options=fate_lmj.var.index,
    description='Fate 2:',
    value=fate_lmj.var.index[1] if len(fate_lmj.var.index) > 1 else fate_lmj.var.index[0]
)
day_dropdown = widgets.Dropdown(
    options=sorted(fate_lmj.obs['day'].unique()),
    description='Day',
    value=2
)

def update_fate_vis(name1, name2, day):
    plt.figure(figsize=(8, 8), facecolor="white")   

    # Extraer las probabilidades de destino para el día seleccionado
    day_mask = fate_lmj.obs['day'] == day
    fate1 = fate_lmj[day_mask, name1].X.flatten()
    fate2 = fate_lmj[day_mask, name2].X.flatten()

    Nrows = len(fate1)
    x = np.zeros(Nrows)
    y = np.zeros(Nrows)
    
    # Matriz con las tres coordenadas del triángulo equilátero
    P = np.array([
        [1.0, 0.0],
        [math.cos(2*math.pi/3), math.sin(2*math.pi/3)],
        [math.cos(4*math.pi/3), math.sin(4*math.pi/3)]
    ])

    # Proyectar el vector tri-dimensional a coordenadas cartesianas (X, Y)
    for i in range(0, Nrows):
        ff = np.array([fate1[i], fate2[i], 1.0 - (fate1[i] + fate2[i])])
        x[i] = (ff @ P)[0]
        y[i] = (ff @ P)[1]

    # Dibujar el fondo del triángulo
    t1 = plt.Polygon(P, color=(0, 0, 0, 0.05), zorder=1)
    plt.gca().add_patch(t1)
    
    # Graficar la nube de puntos celulares
    plt.scatter(x, y, c="#2b8cbe", s=6, alpha=0.6, edgecolors="none", zorder=3)
    plt.scatter(P[:, 0], P[:, 1], c="#e31a1c", s=40, zorder=4)
    
    # Etiquetas de los vértices
    plt.text(P[0, 0] + .05, P[0, 1], name1.replace("_", " ").title(), fontsize=10, fontweight="bold")
    plt.text(P[1, 0] - .10, P[1, 1] + .05, name2.replace("_", " ").title(), fontsize=10, fontweight="bold")
    plt.text(P[2, 0] - .10, P[2, 1] - .10, 'Other Fates', fontsize=10, style="italic")
    
    plt.axis('equal')
    plt.axis('off')
    plt.title(f'Ternary Fate Mapping: {name1.title()} vs. {name2.title()} (Day {day})', fontsize=12, fontweight="bold", pad=15)
    plt.tight_layout()
    plt.show()

# Lanzar los widgets en el cuaderno
widgets.interact(update_fate_vis, name1=fate_dropdown1, name2=fate_dropdown2, day=day_dropdown)

# ============================================================
# Data viewer
# ============================================================

fate_df = fate_lmj.to_df()
fate_df.style.background_gradient(cmap="YlOrRd", vmin=0.0, vmax=1.0).format("{:.3f}")


cell_metadata = pd.read_csv(CELL_METADATA_PATH, index_col=0)

for fate_name in fate_lmj.var_names:

    # --------------------------------------------------------
    # Fate probability
    # --------------------------------------------------------

    fate_values = np.asarray(fate_lmj[:, fate_name].X).flatten()

    # --------------------------------------------------------
    # Cell coordinates
    # --------------------------------------------------------

    plot_data = cell_metadata.loc[fate_lmj.obs_names,["x", "y"]].copy()

    plot_data["fate"] = fate_values

    # --------------------------------------------------------
    # Check
    # --------------------------------------------------------

    print(
        fate_name,
        "n =", len(plot_data),
        "x range =",
        plot_data["x"].min(),
        plot_data["x"].max(),
        "y range =",
        plot_data["y"].min(),
        plot_data["y"].max(),
        "fate range =",
        plot_data["fate"].min(),
        plot_data["fate"].max()
    )

    # --------------------------------------------------------
    # Plot
    # --------------------------------------------------------

    fig, ax = plt.subplots(figsize=(7, 6))

    scatter = ax.scatter(
        plot_data["x"].to_numpy(),
        plot_data["y"].to_numpy(),
        c=plot_data["fate"].to_numpy(),
        s=4,
        cmap="magma",
        vmin=0.1,
        vmax=0.9,
        rasterized=True
    )

    fig.colorbar(
        scatter,
        ax=ax,
        label="Fate probability"
    )

    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")

    ax.set_title(
        f"{fate_name}\n"
        "Fate toward day 9"
    )

    plt.tight_layout()

    plt.show()
    plt.close(fig)

fate_df.to_csv("leishmania_trajectory/results/wot/fate_matrix_table.csv")


