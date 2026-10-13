# =============================================================================
# Extracts the CLL/B-cell compartment from the full PBMC Seurat object,
# exports the coding-gene, zero-filtered expression matrix for pySCENIC,
# and identifies per-timepoint marker genes.
# =============================================================================

source("Scripts/config.R")
source("Scripts/functions.R")

library(Seurat)
library(biomaRt)
library(org.Hs.eg.db)
library(dplyr)
library(scales)

# ── Load and subset ------------------------------------------------------------
load(PBMC_FULL_RDATA)   # loads object `cll_obj` (full PBMC)
cat("Cell types in full PBMC object:\n"); print(table(cll_obj$cell_types))

cll_obj <- subset(cll_obj, subset = cell_types %in% c("CLL", "Bcells"))

cll_obj <- NormalizeData(cll_obj)
cll_obj <- FindVariableFeatures(cll_obj)
cll_obj <- ScaleData(cll_obj)
cll_obj <- RunPCA(cll_obj)
cll_obj <- FindNeighbors(cll_obj, dims = 1:20)
cll_obj <- FindClusters(cll_obj, resolution = 0.5)
cll_obj <- RunUMAP(cll_obj, dims = 1:20)
DimPlot(cll_obj, label = TRUE)

saveRDS(cll_obj, file.path(DATA_DIR, "CLL_only/cll_obj_extracted.rds"))

# ── Export coding-gene, zero-filtered matrix for pySCENIC ---------------------
norm_matrix <- GetAssayData(cll_obj, layer = "data")
norm_matrix <- as.matrix(norm_matrix)

mart      <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_list <- rownames(norm_matrix)

annotations <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = gene_list,
  keytype = "SYMBOL",
  columns = c("SYMBOL", "GENETYPE")
)

coding_genes <- annotations$SYMBOL[annotations$GENETYPE == "protein-coding"]

CLL_normalised_coding      <- norm_matrix[rownames(norm_matrix) %in% coding_genes, ]
CLL_normalised_coding_0Rem <- CLL_normalised_coding[rowSums(CLL_normalised_coding) > 0, ]

write.csv(CLL_normalised_coding_0Rem,
          file.path(DATA_DIR, "CLL_only/CLL_normalized_matrix_coding_0Rem.csv"))


# =============================================================================
# Exploratory checks (timepoint mapping, single-gene expression vs. time)
# =============================================================================

DimPlot(cll_obj, reduction = "umap", group.by = "time_point") +
  ggtitle("Time point mapping") +
  scale_color_manual(values = TIMEPOINT_COLORS)

p1 <- FeaturePlot(cll_obj, "TAF7") + ggtitle("TAF7 expression")
p2 <- FeaturePlot(cll_obj, "RBM3") + ggtitle("RBM3 expression")
p1 | p2

anova(lm(TAF7 ~ time_point, data = FetchData(cll_obj, vars = c("TAF7", "time_point"))))

# ── Per-timepoint marker genes -------------------------------------------------
Idents(cll_obj) <- "time_point"
markers_M0  <- FindMarkers(cll_obj, ident.1 = "M0")
markers_M3  <- FindMarkers(cll_obj, ident.1 = "M3")
markers_M27 <- FindMarkers(cll_obj, ident.1 = "M27")

sig_M0  <- markers_M0  %>% filter(avg_log2FC > 0 & p_val_adj < 0.05) %>% rownames()
sig_M3  <- markers_M3  %>% filter(avg_log2FC > 0 & p_val_adj < 0.05) %>% rownames()
sig_M27 <- markers_M27 %>% filter(avg_log2FC > 0 & p_val_adj < 0.05) %>% rownames()

writeLines(sig_M0,  file.path(RESULTS_DIR, "genes_M0_corrected.txt"))
writeLines(sig_M3,  file.path(RESULTS_DIR, "genes_M3_corrected.txt"))
writeLines(sig_M27, file.path(RESULTS_DIR, "genes_M27_corrected.txt"))

message("Stage 1 complete: CLL compartment extracted, matrix exported for pySCENIC, marker gene lists saved.")
