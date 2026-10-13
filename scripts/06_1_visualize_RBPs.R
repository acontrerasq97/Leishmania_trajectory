library(anndata)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(patchwork)

adata <- read_h5ad("leishmania_trajectory/data/processed/ExprMatrix.preprocessed.h5ad")

df_fates <- read.csv("leishmania_trajectory/results/wot/fate_matrices/fate_matrices_table.csv")
df_expression <- read.csv("leishmania_trajectory/data/processed/ExprMatrix_rbps_haptomonads.csv")
rbps <- read.csv('leishmania_trajectory/results/wot/predictive_rbps/RBPS_haptomonads_9_5.csv') 
rbps_list <- rbps$gene_id

cell_ids <- adata$obs_names
days     <- adata$obs$day


# --- 2. EXTRAER EXPRESIÓN Y METADATOS CELULARES EN MEMORIA ---

df_plot <- data.frame(
  cell_ids = as.character(adata$obs_names),
  days = as.numeric(adata$obs$day),
  expression = as.numeric(as.matrix(adata[, rbps_list]$X))
)


df_plot$weight_haptomonads <- df_fates$haptomonads

# --- 3. REORGANIZAR DATOS EN FORMATO LARGO PARA COMPARAR TRAYECTORIAS ---
# Duplicamos las filas de las células para asociarlas a cada trayectoria con su respectivo peso
df_trajectories <- df_plot %>%
  pivot_longer(
    cols = c(weight_haptomonads),
    names_to = "Trajectory",
    names_prefix = "weight_",
    values_to = "Probability"
  ) %>%
  # Añadimos un sutil ruido aleatorio (jitter) en el eje X (day) para dispersar las células
  # y que no queden todas aplastadas en columnas estrictas (0, 5, 10, 15)
  mutate(pseudotime = days + runif(n(), -0.5, 0.5))

# --- 4. GRAFICAR LAS CELULAS Y EL AJUSTE SUAVIZADO GAM ---
trajectory_plot <- ggplot(df_trajectories, aes(x = pseudotime, y = expression, color = Trajectory)) +
  # 1. Graficar células individuales. El tamaño y opacidad dependen de su probabilidad en la ruta
  geom_point(aes(alpha = Probability, size = Probability)) +
  
  # 2. Ajustar el modelo suavizado continuo (GAM Spline) sobre las células para trazar la tendencia
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), linewidth = 1.5, se = TRUE) +
  
  # 3. Configuración estética y de colores para publicación
  scale_color_manual(values = c("haptomonads" = "#E63946", "procyclics" = "#457B9D")) +
  scale_alpha_continuous(range = c(0.05, 0.6), guide = "none") +
  scale_size_continuous(range = c(0.5, 2.0), guide = "none") +
  
  labs(
    title = paste("Single-Cell Kinetics of", rbps_list, "Across WOT Trajectories"),
    x = "Developmental Timeline (Days with Jitter)",
    y = "Normalized Expression Level",
    color = "Lineage Trajectory"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", h = 0.5),
    legend.position = "top",
    legend.title = element_text(face = "bold")
  )

# Desplegar el gráfico definitivo
print(trajectory_plot)

