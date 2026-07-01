# Load the libraries
library(Seurat)
library(tidyverse)

# Load the regulon matrix
auc_mat <- read.csv("Results/CLL_only/AUcell_CLL_only.csv", header = T, row.names = 1, sep = "\t")

# Remove the ... in the column names
colnames(auc_mat) <- gsub("\\.\\.\\.", "", colnames(auc_mat))
cell_names <- rownames(auc_mat)

# Load the seurat object and extract the gene expression matrix from there 
seurat_obj <- LoadSeuratRds("Data/CLL_only/cll_obj_tfAct_AUCell.rds")
gene_expr <- GetAssayData(
  seurat_obj,
  assay = "RNA",
  layer  = "data"
)
gene_expr <- as.matrix(gene_expr) %>% t() %>% as.data.frame()

# Make sure we have the same cells in the two objects 

common_cells <- intersect(rownames(gene_expr), rownames(auc_mat))
gene_expr <- gene_expr[common_cells, ]
auc_mat   <- auc_mat[common_cells, ]

# Pseudobulking for each time point 
timepoints <- c("M0", "M3", "M27")

expr_pb <- sapply(timepoints, function(tp) {
  cells <- rownames(seurat_obj@meta.data)[seurat_obj$time_point == tp]
  colMeans(gene_expr[cells, , drop = FALSE])
})

auc_pb <- sapply(timepoints, function(tp) {
  cells <- rownames(seurat_obj@meta.data)[seurat_obj$time_point == tp]
  colMeans(auc_mat[cells, , drop = FALSE])
})

# Load the GRN 
net_csv <- read.csv("Results/CLL_only/SCENIC_GRN_edges_allfilters.csv", header = T)

# From the GRN select the Source from the AUCell matrix and the target from the gene expression matrix
grn_edges <- net_csv[, c("Source", "Target", "Weight")]
grn_edges <- grn_edges[
  grn_edges$Source %in% rownames(auc_pb) &
    grn_edges$Target %in% rownames(expr_pb),
]

# Calculate the Spearman correlation 
grn_edges$rho <- mapply(
  function(tf, gene) {
    cor(
      auc_pb[tf, ],
      expr_pb[gene, ],
      method = "spearman"
    )
  },
  grn_edges$Source,
  grn_edges$Target
)

# Assign correlation values to edge signs 
grn_edges$sign <- with(grn_edges, ifelse(
  rho >=  0.3, "activation",
  ifelse(rho <= -0.3, "inhibition", "activation")
))

# Save the signed GRN 
write.csv(grn_edges, "Results/CLL_only/SCENIC_GRN_edges_allfilters_signed.csv", quote = FALSE)
