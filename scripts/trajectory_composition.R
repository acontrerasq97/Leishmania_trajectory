library(anndata)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# --- 1. LOAD YOUR EXACT INPUT FILES ---
adata    <- read_h5ad("leishmania_trajectory/data/processed/ExprMatrix.preprocessed.h5ad")
df_fates <- read.csv("leishmania_trajectory/results/wot/fate_matrices/fate_matrices_table.csv")

gmt_lines <- readLines("leishmania_trajectory/data/metadata/cell_sets.gmt")

#Calculate the entropy per cell: low values ​​imply a single, consolidated fates
df_fates$entropy <- apply(df_fates[, c("haptomonads", "late_metacyclics", "leptomonads_2")], 1, function(p) {
  p <- p[p > 0]
  -sum(p * log2(p))
})

#Indice de Shannon
df_meta$entropy <- df_fates$entropy
df_entropy_summary <- df_meta %>%
  group_by(cell_type) %>%
  summarise(
    mean_entropy = mean(entropy, na.rm = TRUE),
    sd_entropy   = sd(entropy, na.rm = TRUE),
    total_cells  = n()
  ) %>%
  arrange(mean_entropy)
df_entropy_summary


cell_to_stage_list <- list()

for (line in gmt_lines) {
  # Split the tab-separated line
  tokens <- unlist(strsplit(line, "\t"))
  if (length(tokens) > 2) {
    stage_name <- tokens[1]
    # The remaining tokens (from index 3 onwards) are the cell IDs
    cell_ids_in_stage <- tokens[3:length(tokens)]
    # Trim whitespaces just in case
    cell_ids_in_stage <- trimws(cell_ids_in_stage)
    # Exclude any trailing empty strings
    cell_ids_in_stage <- cell_ids_in_stage[cell_ids_in_stage != ""]
    
    # Store in a temporary named list
    cell_to_stage_list[[stage_name]] <- cell_ids_in_stage
  }
}

df_gmt_lookup <- stack(cell_to_stage_list)
colnames(df_gmt_lookup) <- c("cell_id", "cell_type")


df_meta_base <- data.frame(
  cell_id   = as.character(adata$obs_names$to_numpy()),
  day       = as.factor(as.numeric(as.character(adata$obs$day)))
)

df_meta <- df_meta_base %>% left_join(df_gmt_lookup, by = "cell_id")


# Append the pre-calculated WOT probabilities from your fates matrix
df_meta$weight_haptomonads    <- df_fates$haptomonads
df_meta$weight_metacyclics   <- df_fates$late_metacyclics 
df_meta$weight_leptomonads_2   <- df_fates$leptomonads_2
df_meta$weight_procyclics    <- df_fates$procyclics
df_meta$weight_early_metacyclics   <- df_fates$early_metacyclics 
df_meta$weight_leptomonads_1   <- df_fates$leptomonads_1
df_meta$weight_procyclics_to_nectomonads    <- df_fates$procyclics_to_nectomonads
df_meta$weight_nectomonads   <- df_fates$nectomonads 
df_meta$weight_leptomonads_0   <- df_fates$leptomonads_0
df_meta$weight_amastigotes   <- df_fates$amastigotes


# --- 3. PROPORTIONAL COMPOSITION FUNCTION ---
get_trajectory_composition <- function(meta_df, weight_column, lineage_label) {
  composition <- meta_df %>%
    # Use the specific trajectory weight for the calculations
    mutate(w = .data[[weight_column]]) %>%
    # Group by time point and biological stage
    group_by(day, cell_type) %>%
    # Sum the fate probabilities to get the weighted total mass of cells
    summarise(weighted_cell_count = sum(w), .groups = 'drop') %>%
    mutate(Lineage = lineage_label)
  
  return(composition)
}

# --- 4. PROCESS BOTH LINEAGES ---
comp_haptomonads  <- get_trajectory_composition(df_meta, "weight_haptomonads", "haptomonads Trajectory")
comp_metacyclics  <- get_trajectory_composition(df_meta, "weight_metacyclics", "late-metacyclics Trajectory")
comp_leptomonads2 <- get_trajectory_composition(df_meta, "weight_leptomonads_2", "leptomonads_2 Trajectory")
comp_leptomonads0 <- get_trajectory_composition(df_meta, "weight_leptomonads_0", "leptomonads_0 Trajectory")
comp_leptomonads1 <- get_trajectory_composition(df_meta, "weight_leptomonads_1", "leptomonads_1 Trajectory")
comp_nectopmonads<- get_trajectory_composition(df_meta, "weight_nectomonads", "nectomonads Trajectory")

# --- 5. GENERATE STACKED BAR PLOTS (MATCHING YOUR IMAGE) ---
# Custom plotting theme to match the clean publication style from your image
theme_publication_bar <- function() {
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", h = 0.5, margin = margin(b = 15)),
    axis.text.x = element_text(angle = 30, hj = 1),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )
}

my_colors <- c(
  # --- Estadios Tempranos e Intermedios Base ---
  "amastigotes"               = "#56B4E9",  # Azul claro (Excelente visibilidad inicial, alto contraste)
  "procyclics"                = "#0072B2",  # Azul oscuro (Mismo tono, pero caída drástica en luminosidad)
  "procyclics_to_nectomonads" = "#009E73",  # Verde azulado/esmeralda (Apto para daltónicos, evita el verde pasto)
  
  # --- Serie de Leptomónadas (Gradiente de diferenciación secuencial) ---
  "leptomonads_0"             = "#351A44",  # Morado berenjena muy oscuro (Aporta contraste por sombra)
  "leptomonads_1"             = "#785EF0",  # Violeta medio (Paleta accesible de IBM, alta discriminación)
  "leptomonads_2"             = "#a30ed4",  # Lavanda claro/pastel (Contraste por brillo con los anteriores)
  
  # --- Estadios Avanzados y de Transición de Lípidos/Estructura ---
  "nectomonads"               = "#E69F00",  # Naranja ocre brillante (Paleta Okabe-Ito, muy distinguible de azules)
  
  # --- Destino Terminal: Metacíclicos (Camino Alternativo) ---
  "early_metacyclics"         = "#F0E442",  # Amarillo suave (Destaca inmediatamente por alta luminosidad)
  "late_metacyclics"          = "#CC79A7",  # Rosa pálido/magento (Diferenciación total con la ruta haptomónada)
  
  # --- Destino Terminal: Haptomónadas (Tu ruta de interés) ---
  "haptomonads"               = "#D55E00"  # Rojo bermellón / Naranja rojizo (Fuerte impacto visual al final de la barra)
)

# Plot 1: Haptomonads Composition over Time
p_hapto <- ggplot(comp_leptomonads2, aes(x = day, y = weighted_cell_count, fill = cell_type)) +
  geom_bar(stat = "identity", position = "fill", color = "white", linewidth = 0.3) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  scale_fill_manual(values = my_colors, name = "Cell Stage") +
  labs(title = "Lineage: Haptomonads", x = "Developmental Timeline (days)", y = "Frequency") +
  theme_publication_bar()

# Plot 2: Late-Metacyclics Composition over Time
p_meta <- ggplot(comp_metacyclics, aes(x = day, y = weighted_cell_count, fill = cell_type)) +
  geom_bar(stat = "identity", position = "fill", color = "white", linewidth = 0.3) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  scale_fill_manual(values = my_colors, name = "Cell Stage") +
  labs(title = "Lineage: Late-Metacyclics", x = "Developmental Timeline (days)", y = "Frequency") +
  theme_publication_bar() +
  theme(axis.title.y = element_blank()) # Hide Y axis label on the second plot for clarity

# Combine both plots side-by-side with a single shared legend using patchwork
final_composition_plot <- (p_hapto + p_meta) + plot_layout(guides = 'collect') & theme(legend.position = "right")

# Save the high-resolution figure
ggsave("leishmania_trajectory/figures/wot/wot_trajectories_stacked_composition.png", plot = final_composition_plot, width = 12, height = 7, dpi = 300)
print(final_composition_plot)


# --- FUNCIÓN PARA GENERAR EL GRÁFICO DESGLOSADO POR ESTADIO ---
generate_facet_plot <- function(composition_df, lineage_title) {
  ggplot(composition_df, aes(x = day, y = weighted_cell_count, fill = cell_type)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3, show.legend = FALSE) +
    
      scale_y_continuous(
      trans = "log10", 
      labels = scales::label_scientific(digits = 1), 
      name = "Weighted Cell Mass (Log10 Scale)"
    ) +
    
    scale_fill_manual(values = my_colors) +
    facet_wrap(~cell_type, scales = "free_y", ncol = 2) + 
    labs(title = paste("Lineage:", lineage_title), x = "Developmental Timeline") +
    
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", h = 0.5, size = 14, margin = margin(b = 10)),
      strip.text = element_text(face = "bold", size = 10),
      axis.text.x = element_text(angle = 0, h = 0.5),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "#F2F2F2")
    )
}

p_final_hapto <- generate_facet_plot(comp_haptomonads, "Haptomonads")
p_final_meta  <- generate_facet_plot(comp_metacyclics, "Late-Metacyclics")

final_face_plot <- (p_final_hapto + p_final_meta) + plot_layout(guides = 'collect') 

ggsave("leishmania_trajectory/figures/wot/wot_composition_facets.png", plot = final_face_plot, width = 9, height = 11, dpi = 300)

print(final_face_plot)

