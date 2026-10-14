library(Seurat)
library(plotly)
library(htmlwidgets)
library(viridisLite) # Para asegurar el acceso a paletas científicas
library(plotly)
library(scales)

# Assumes `seu` is a Seurat object that already went through
# NormalizeData, FindVariableFeatures, ScaleData, RunPCA and FindClusters.

seu <- readRDS("~/Documents/LmjRy/LmjRyDevelopment.rds")

seu <- FindSubCluster(
  object = seu,
  cluster = "leptomonads",              # El cluster exacto que quieres dividir
  graph.name = "SCT_snn",              # El nombre de tu grafo de vecinos de SCTransform
  subcluster.name = "subclusters_leptomonads", # Nueva columna que se creará en meta.data
  resolution = 0.1                
)

# 1. Compute a 3-component UMAP stored under its own reduction name
seu <- RunUMAP(
  seu,
  dims           = 1:30,
  n.components   = 3L,
  reduction.name = "umap3d",
  reduction.key  = "UMAP3D_",
  seed.use       = 42
)

plot_umap_3d <- function(
                        seu, 
                        color_by = "subclusters_leptomonads", 
                        reduction = "umap3d", 
                        point_size = 1, 
                        opacity = 0.85, 
                        palette = NULL) {
  
  emb <- as.data.frame(Embeddings(seu, reduction = reduction))

  colnames(emb) <- c("UMAP_1", "UMAP_2", "UMAP_3")
  
  if (color_by %in% colnames(seu@meta.data)) {
    emb[["color"]] <- seu@meta.data[[color_by]]
  } else if (color_by %in% rownames(seu)) {
    emb[["color"]] <- FetchData(seu, vars = color_by)[, 1]
  } else {
    stop("`color_by` was not found in metadata or features.")
  }
  emb[["cell"]] <- rownames(emb)
  is_discrete <- is.factor(emb[["color"]]) || is.character(emb[["color"]])
  if (is_discrete) {
    emb[["color"]] <- factor(emb[["color"]])
    if (is.null(palette)) {
      okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000", "#999999")      
      if (nlevels(emb[["color"]]) <= length(okabe_ito)) {
        palette <- okabe_ito[1:nlevels(emb[["color"]])]
      } else {
        palette <- viridisLite::turbo(nlevels(emb[["color"]]))
      }
    }
    marker <- list(size = point_size, opacity = opacity)
    colors_arg <- palette
  } else {
    marker <- list(size = point_size, 
      opacity = opacity, 
      color = emb[["color"]], 
      colorscale = "Cividis", 
      showscale = TRUE, 
      colorbar = list(title = list(text = color_by, 
        font = list(size = 12)), 
        tickfont = list(size = 10), 
        thickness = 15))
    colors_arg <- NULL
  }
  p <- if (is_discrete) {
    plot_ly(data = emb, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
      color = ~color, 
      colors = colors_arg, 
      type = "scatter3d", 
      mode = "markers", 
      marker = marker, 
      text = ~cell, 
      hoverinfo = "text+name")
  } else {
    plot_ly(data = emb, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
      type = "scatter3d", 
      mode = "markers", 
      marker = marker, 
      text = ~paste0("Célula: ", cell, "<br>", color_by, ": ", round(color, 2)), 
      hoverinfo = "text")
  }
  p %>% 
    layout(title = list(text = paste("UMAP 3D:", color_by), 
  font = list(size = 14, 
    face = "bold")), 
    scene = list(xaxis = list(title = "UMAP 1", 
    backgroundcolor = "#f0f0f0", 
    showbackground = TRUE, 
    gridcolor = "#ffffff", 
    zerolinecolor = "#ffffff"), 
    yaxis = list(title = "UMAP 2", 
    backgroundcolor = "#f0f0f0", 
    showbackground = TRUE, 
    gridcolor = "#ffffff", 
    zerolinecolor = "#ffffff"), 
    zaxis = list(title = "UMAP 3", 
    backgroundcolor = "#e5e5e5", 
    showbackground = TRUE, 
    gridcolor = "#ffffff", 
    zerolinecolor = "#ffffff"), 
    camera = list(eye = list(x = 1.6, y = 1.6, z = 1.3))), 
    legend = list(itemsizing = "constant", font = list(size = 11)), 
    margin = list(l = 0, r = 0, b = 0, t = 40))
}




# 3. Examples
p_clusters <- plot_umap_3d(seu, color_by = "subclusters_leptomonads")
p_gene     <- plot_umap_3d(seu, color_by = "CD3E")

p_clusters

# 4. Save as a standalone interactive HTML file
saveWidget(p_clusters, "leishmania_trajectory/figures/wot/umap_3d_clusters.html", selfcontained = TRUE)
