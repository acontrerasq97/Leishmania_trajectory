library(Seurat)
library(plotly)
library(htmlwidgets)

# Assumes `seu` is a Seurat object that already went through
# NormalizeData, FindVariableFeatures, ScaleData, RunPCA and FindClusters.

seu <- readRDS("~/Documents/LmjRy/LmjRyDevelopment.rds")

seu <- FindSubCluster(
  object = seu,
  cluster = "leptomonads",              # El cluster exacto que quieres dividir
  graph.name = "SCT_snn",              # El nombre de tu grafo de vecinos de SCTransform
  subcluster.name = "subclusters_leptomonads", # Nueva columna que se creará en meta.data
  resolution = 0.2                
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

# 2. Interactive 3D plot colored by any metadata column or gene
plot_umap_3d <- function(seu,
                         color_by  = "clusters",
                         reduction = "umap3d",
                         point_size = 2,
                         opacity    = 0.8,
                         palette    = NULL) {

  emb <- as.data.frame(Embeddings(seu, reduction = reduction))
  colnames(emb) <- c("UMAP_1", "UMAP_2", "UMAP_3")

  if (color_by %in% colnames(seu@meta.data)) {
    emb$color <- seu@meta.data[[color_by]]
  } else if (color_by %in% rownames(seu)) {
    emb$color <- FetchData(seu, vars = color_by)[, 1]
  } else {
    stop("`color_by` was not found in metadata or features.")
  }
  emb$cell <- rownames(emb)

  is_discrete <- is.factor(emb$color) || is.character(emb$color)

  if (is_discrete) {
    emb$color <- factor(emb$color)
    if (is.null(palette)) {
      palette <- scales::hue_pal()(nlevels(emb$color))
    }
    marker <- list(size = point_size, opacity = opacity)
    colors_arg <- palette
  } else {
    marker <- list(
      size = point_size, opacity = opacity,
      color = emb$color, colorscale = "Viridis", showscale = TRUE,
      colorbar = list(title = color_by)
    )
    colors_arg <- NULL
  }

  p <- if (is_discrete) {
    plot_ly(
      data = emb,
      x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,
      color = ~color, colors = colors_arg,
      type = "scatter3d", mode = "markers",
      marker = marker,
      text = ~cell, hoverinfo = "text+name"
    )
  } else {
    plot_ly(
      data = emb,
      x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,
      type = "scatter3d", mode = "markers",
      marker = marker,
      text = ~paste0(cell, "<br>", color_by, ": ", round(color, 2)),
      hoverinfo = "text"
    )
  }

  p %>% layout(
    title = color_by,
    scene = list(
      xaxis = list(title = "UMAP 1"),
      yaxis = list(title = "UMAP 2"),
      zaxis = list(title = "UMAP 3")
    ),
    legend = list(itemsizing = "constant")
  )
}

# 3. Examples
p_clusters <- plot_umap_3d(seu, color_by = "clusters")
p_gene     <- plot_umap_3d(seu, color_by = "CD3E")

p_clusters

# 4. Save as a standalone interactive HTML file
saveWidget(p_clusters, "leishmania_trajectory/figures/wot/umap_3d_clusters.html", selfcontained = TRUE)
