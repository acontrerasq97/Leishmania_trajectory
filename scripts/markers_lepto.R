library(Seurat)
library(dplyr)
library(ggplot2)
library(zellkonverter)

sce <- readH5AD("data/processed/leptomonads_subset.h5ad")

seurat_obj <- as.Seurat(sce, counts = "X", data = "X")
Idents(seurat_obj) <- "original_subpopulation"
table(Idents(seurat_obj))
# 1. Filtrar solo las células correspondientes a leptomonads
# Asumiendo que la columna de metadatos se llama 'cell_subpopulation'
lepto_cells <- subset(seurat_obj, idents = c("leptomonads_0", "leptomonads_1", "leptomonads_2"))

# Verificar el número de células por grupo
table(Idents(lepto_cells))

# 2. Encontrar marcadores específicos (up-regulated) para cada subgrupo
lepto_markers <- FindAllMarkers(
  lepto_cells,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

# Guardar la tabla completa para tener el respaldo numérico
write.csv(lepto_markers, "./leptomonad_subcluster_markers_Seurat.csv", row.names = FALSE)

# Ver cuántos marcadores pasa cada grupo con un p-adj < 0.05
lepto_markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  summarise(n_genes = n())
