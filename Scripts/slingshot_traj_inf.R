library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(scater)
library(ggplot2)

seu <- LoadSeuratRds("Data/CLL_only/cll_obj_tfAct_AUCell.rds")
DefaultAssay(seu) <- "RNA"
sce <- as.SingleCellExperiment(seu)

reducedDim(sce, "UMAP") <- seu@reductions$umap@cell.embeddings
umap <- reducedDim(sce, "UMAP")
umap[,2] <- -umap[,2]  # flip vertically
reducedDim(sce, "UMAP") <- umap

sce$cluster <- factor(seu$seurat_clusters)

table(sce$cluster, sce$time_point)

start_cluster <- "3"

sce <- slingshot(
  sce,
  clusterLabels = "cluster",
  reducedDim = "UMAP",
  start.clus = start_cluster
)

plot(
  reducedDim(sce, "UMAP"),
  col = as.numeric(sce$cluster),
  pch = 16,
  asp = 1
)

lines(
  SlingshotDataSet(sce),
  lwd = 2,
  col = "black"
)

time_cols <- c(
  "M0"  = "#66C2A5",
  "M3"  = "#FC8D62",
  "M27" = "#8DA0CB"
)

plot(
  reducedDim(sce, "UMAP"),
  col = time_cols[sce$time_point],
  pch = 16,
  asp = 1,
  main = "Slingshot trajectories"
)

# Add Slingshot trajectories
lines(SlingshotDataSet(sce), lwd = 2, col = "royalblue")

# Add legend
legend(
  "topright",
  legend = names(time_cols),
  col = time_cols,
  pch = 16,
  bty = "n"
)
