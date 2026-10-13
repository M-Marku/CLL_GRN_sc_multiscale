# =============================================================================
# CopyKAT CNV inference (T-cell reference), optimal-k diagnostics, clone
# identification, clonal tree, resistance-clone markers, ancestral state
# reconstruction. Produces manuscript Figure 4 (a-e).
# =============================================================================

source("config.R")
source("functions.R")

library(Seurat)
library(tidyverse)
library(ape)
library(phytools)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(cluster)
library(copykat)

# Local overrides (this stage uses its own output dir, kept separate from the
# other stages' RESULTS_DIR/PLOTS_DIR since Fig 4 outputs are self-contained)
OUT_DIR <- PHYLO_OUT_DIR
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

N_HVG    <- 2000
N_PCS    <- 20
K_RANGE  <- 2:8
K_CLONES <- 3   # revisit after Section 2 diagnostics


# =============================================================================
# SECTION 1 — CNV INFERENCE WITH COPYKAT (T-CELL REFERENCE)              
# =============================================================================

message("\n=== SECTION 1: CNV inference (T-cell reference) ===\n")

seu <- readRDS(CLL_ONLY_RDS)
DefaultAssay(seu) <- "RNA"
cat("CLL object: ", ncol(seu), "cells,", nrow(seu), "genes\n")
cat("Timepoints:\n"); print(table(seu[[TIMEPOINT_COL]]))

load(PBMC_FULL_RDATA)   # loads object `cll_obj` (full PBMC, used only for T-cell reference here)
DefaultAssay(cll_obj) <- "RNA"
stopifnot(CELLTYPE_COL %in% colnames(cll_obj@meta.data))
cat("\nCell types available in the full PBMC object:\n")
print(table(cll_obj[[CELLTYPE_COL]][, 1]))

ref_cells <- colnames(cll_obj)[cll_obj[[CELLTYPE_COL]][, 1] %in% REF_CELLTYPES]
cat("\nT-cell reference cells found:", length(ref_cells), "\n")
stopifnot(length(ref_cells) > 20)
stopifnot(length(intersect(ref_cells, colnames(seu))) == 0)

counts_cll <- GetAssayData(seu,     assay = "RNA", layer = "counts")
counts_ref <- GetAssayData(cll_obj, assay = "RNA", layer = "counts")[, ref_cells]

shared_genes <- intersect(rownames(counts_cll), rownames(counts_ref))
cat("Shared genes between CLL and reference matrices:", length(shared_genes), "\n")

raw_counts <- cbind(counts_cll[shared_genes, ], counts_ref[shared_genes, ])
cat("Combined matrix:", ncol(raw_counts), "cells (",
    ncol(counts_cll), "CLL candidates +", length(ref_cells), "T-cell reference )\n")

RUN_COPYKAT <- TRUE

if (RUN_COPYKAT) {
  ck <- copykat(
    rawmat          = as.matrix(raw_counts),
    id.type         = "S",
    ngene.chr       = 5,
    win.size        = 25,
    KS.cut          = 0.1,
    LOW.DR          = 0.02,
    UP.DR           = 0.02,
    sam.name        = "CLL_ibrutinib_Tref",
    distance        = "euclidean",
    norm.cell.names = ref_cells,
    output.seg      = TRUE,
    plot.genes      = FALSE,
    n.cores         = 4
  )
  saveRDS(ck, file.path(OUT_DIR, "copykat_result_Tref.rds"))
} else {
  ck <- readRDS(file.path(OUT_DIR, "copykat_result_Tref.rds"))
}

cnv_mat <- t(ck$CNAmat[, -c(1:3)])
rownames(cnv_mat) <- gsub("\\.1$", "-1", rownames(cnv_mat))

pred_df <- ck$prediction
rownames(pred_df) <- pred_df$cell.names

shared_cells <- Reduce(intersect, list(colnames(seu), rownames(pred_df), rownames(cnv_mat)))
cnv_mat <- cnv_mat[shared_cells, ]
pred_df <- pred_df[shared_cells, ]

seu$copykat_pred <- pred_df[colnames(seu), "copykat.pred"]

cat("\nCopyKAT classification of CLL cells (T-cell reference):\n")
print(table(seu$copykat_pred, useNA = "ifany"))

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, nfeatures = N_HVG)
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = N_PCS)
seu <- RunUMAP(seu, dims = 1:N_PCS)

p_umap_ck <- DimPlot(seu, group.by = "copykat_pred") +
  ggtitle("CopyKAT classification (T-cell reference)") +
  scale_color_manual(values = c("aneuploid" = "#e74c3c", "diploid" = "royalblue"))
ggsave(file.path(OUT_DIR, "Fig4a_umap_copykat_pred.svg"), p_umap_ck, width = 6, height = 5)

saveRDS(seu, file.path(OUT_DIR, "seu_with_copykat_pred.rds"))
message("Section 1 complete.")


# =============================================================================
# SECTION 2 — OPTIMAL NUMBER OF CLONES                                   
# =============================================================================

message("\n=== SECTION 2: Optimal number of clones ===\n")

tumour_cells <- rownames(pred_df)[pred_df$copykat.pred == "aneuploid"]
cnv_tumour   <- cnv_mat[tumour_cells, ]
cat("Aneuploid (tumour) cells entering clustering:", length(tumour_cells), "\n")

cnv_dist <- dist(cnv_tumour, method = "euclidean")
saveRDS(cnv_dist, file.path(OUT_DIR, "cnv_dist.rds"))   # cache — O(n^2), see notes

hc_cnv <- hclust(cnv_dist, method = "ward.D2")

## 2.1 Elbow (centroid-based WSS — fast, avoids repeated O(n^2) dist() calls)
wss_by_k <- sapply(K_RANGE, function(k) {
  cl <- cutree(hc_cnv, k = k)
  sum(sapply(unique(cl), function(cli) {
    idx <- which(cl == cli)
    if (length(idx) < 2) return(0)
    centroid <- colMeans(cnv_tumour[idx, , drop = FALSE])
    sum(scale(cnv_tumour[idx, , drop = FALSE], center = centroid, scale = FALSE)^2)
  }))
})
names(wss_by_k) <- K_RANGE
saveRDS(wss_by_k, file.path(OUT_DIR, "wss_by_k.rds"))

p_elbow <- ggplot(data.frame(k = K_RANGE, wss = wss_by_k), aes(k, wss)) +
  geom_line() + geom_point(size = 2) +
  scale_x_continuous(breaks = K_RANGE) +
  labs(title = "Elbow plot", x = "Number of clones (k)", y = "Within-cluster SS") +
  theme_classic(base_size = 13)
ggsave(file.path(OUT_DIR, "S_optimal_k_elbow.svg"), p_elbow, width = 5, height = 4)

## 2.2 Average silhouette width
sil_by_k <- sapply(K_RANGE, function(k) {
  cl <- cutree(hc_cnv, k = k)
  mean(silhouette(cl, cnv_dist)[, "sil_width"])
})
p_sil <- ggplot(data.frame(k = K_RANGE, sil = sil_by_k), aes(k, sil)) +
  geom_line() + geom_point(size = 2) +
  scale_x_continuous(breaks = K_RANGE) +
  labs(title = "Average silhouette width", x = "Number of clones (k)",
       y = "Mean silhouette width") +
  theme_classic(base_size = 13)
ggsave(file.path(OUT_DIR, "S_optimal_k_silhouette.svg"), p_sil, width = 5, height = 4)
best_k_sil <- K_RANGE[which.max(sil_by_k)]

## 2.3 Gap statistic — subsampled for tractability (see chat notes: full-data
## version can run overnight on ~9k cells). Drop the subsampling block if you
## want the full-data version and have time to let it run.
set.seed(42)
N_SUB <- min(2000, nrow(cnv_tumour))
sub_idx <- sample(nrow(cnv_tumour), N_SUB)
cnv_sub <- cnv_tumour[sub_idx, ]
cat("Gap statistic computed on a random subsample of", N_SUB,
    "of", nrow(cnv_tumour), "aneuploid cells (for tractability)\n")

FUNcluster_ward <- function(x, k) {
  list(cluster = cutree(hclust(dist(x, method = "euclidean"), method = "ward.D2"), k = k))
}
gap_res <- clusGap(cnv_sub, FUNcluster = FUNcluster_ward, K.max = max(K_RANGE), B = 20)
saveRDS(gap_res, file.path(OUT_DIR, "gap_stat.rds"))

gap_tab <- as.data.frame(gap_res$Tab)
gap_tab$k <- as.integer(rownames(gap_tab))

p_gap <- ggplot(gap_tab, aes(k, gap)) +
  geom_errorbar(aes(ymin = gap - SE.sim, ymax = gap + SE.sim), width = 0.15) +
  geom_line() + geom_point(size = 2) +
  scale_x_continuous(breaks = gap_tab$k) +
  labs(title = "Gap statistic (subsampled, B = 20 bootstraps)",
       x = "Number of clones (k)", y = "Gap") +
  theme_classic(base_size = 13)
ggsave(file.path(OUT_DIR, "S_optimal_k_gapstat.svg"), p_gap, width = 5, height = 4)

best_k_gap <- maxSE(gap_tab$gap, gap_tab$SE.sim, method = "firstSEmax")

cat("\n--- Optimal-k diagnostics summary ---\n")
cat("Best k (silhouette, max width)      :", best_k_sil, "\n")
cat("Best k (gap statistic, 1-SE rule)   :", best_k_gap, "\n")
cat("k used downstream (K_CLONES)        :", K_CLONES, "\n")
if (best_k_sil != K_CLONES || best_k_gap != K_CLONES) {
  warning("Diagnostics do not unanimously agree with K_CLONES = ", K_CLONES,
          ". Inspect the elbow/silhouette/gap plots before finalising k for the paper.")
}
message("Section 2 complete.")


# =============================================================================
# SECTION 3 — CLONE IDENTIFICATION & PROPORTIONS ACROSS TIMEPOINTS   
# =============================================================================

message("\n=== SECTION 3: Clone identification ===\n")

clone_labels <- cutree(hc_cnv, k = K_CLONES)
names(clone_labels) <- tumour_cells

seu$clone <- NA
seu$clone[names(clone_labels)] <- paste0("Clone_", clone_labels)
seu$clone[seu$copykat_pred == "diploid"] <- "Diploid_ref"

p_umap_clone <- DimPlot(seu, group.by = "clone", label = TRUE, label.size = 5) +
  ggtitle("Clone identity") +
  scale_colour_manual(values = CLONE_COLORS, na.value = "#BDC3C7")
ggsave(file.path(OUT_DIR, "Fig4b_umap_clones.svg"), p_umap_clone, width = 6, height = 5)

clone_time_tbl <- seu@meta.data %>%
  filter(!is.na(clone), clone != "Diploid_ref") %>%
  count(!!sym(TIMEPOINT_COL), clone) %>%
  group_by(!!sym(TIMEPOINT_COL)) %>%
  mutate(prop = n / sum(n))

p_bar <- ggplot(clone_time_tbl, aes(x = .data[[TIMEPOINT_COL]], y = prop, fill = clone)) +
  geom_col() +
  scale_fill_manual(values = CLONE_COLORS[c("Clone_1", "Clone_2", "Clone_3")]) +
  labs(title = "Clone proportions across timepoints", x = "Timepoint", y = "Proportion") +
  theme_classic() +
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12),
        legend.text = element_text(size = 12), plot.title = element_text(size = 14, face = "bold"),
        legend.title = element_blank(), legend.position = "top")
ggsave(file.path(OUT_DIR, "Fig4c_clone_proportions.svg"), p_bar, width = 5, height = 5)

clone_prop_wide <- clone_time_tbl %>%
  ungroup() %>% dplyr::select(-n) %>%
  pivot_wider(names_from = clone, values_from = prop, values_fill = 0) %>%
  column_to_rownames(TIMEPOINT_COL)

cat("\nClone proportions by timepoint:\n"); print(round(clone_prop_wide, 3))

saveRDS(seu, file.path(OUT_DIR, "seu_with_clones.rds"))
message("Section 3 complete.")


# =============================================================================
# SECTION 4 — CLONAL TREE (NJ, BOOTSTRAP)                               
# =============================================================================

message("\n=== SECTION 4: Clonal tree ===\n")

clone_ids   <- paste0("Clone_", 1:K_CLONES)
clone_cents <- sapply(clone_ids, function(cl) {
  cells_cl <- names(clone_labels)[clone_labels == as.integer(sub("Clone_", "", cl))]
  colMeans(cnv_mat[cells_cl, , drop = FALSE])
})
clone_cents <- t(clone_cents)

germline  <- matrix(0, nrow = 1, ncol = ncol(clone_cents),
                     dimnames = list("Germline", colnames(clone_cents)))
clone_all <- rbind(germline, clone_cents)
clone_dist <- dist(clone_all, method = "euclidean")

tree_nj <- nj(clone_dist)
tree_nj <- root(tree_nj, outgroup = "Germline", resolve.root = TRUE)

n_boot <- 100
boot_trees <- vector("list", n_boot)
for (i in seq_len(n_boot)) {
  boot_cols       <- sample(ncol(clone_all), replace = TRUE)
  boot_dist       <- dist(clone_all[, boot_cols], method = "euclidean")
  boot_trees[[i]] <- nj(boot_dist)
}
class(boot_trees) <- "multiPhylo"
bs_props <- prop.clades(tree_nj, boot_trees)
tree_nj$node.label <- paste0(round(bs_props / n_boot * 100), "%")

svg(file.path(OUT_DIR, "Fig4d_tree_NJ_bootstrap.svg"), width = 6, height = 5)
plot(tree_nj, show.node.label = TRUE,
     tip.color = c("grey40", brewer.pal(K_CLONES, "Set2")),
     main = "NJ tree, rooted (Germline), 100 bootstraps")
add.scale.bar()
dev.off()

message("Section 4 complete.")


# =============================================================================
# SECTION 5 — TRANSCRIPTOMIC CHARACTERISATION OF THE RESISTANCE CLONE   
# =============================================================================

message("\n=== SECTION 5: Clone markers & volcano plot ===\n")

Idents(seu) <- "clone"
seu_tumour <- subset(seu, subset = clone %in% clone_ids)

markers_list <- lapply(clone_ids, function(cl) {
  markers <- FindMarkers(seu_tumour, ident.1 = cl,
                          min.pct = 0.1, logfc.threshold = 0.25, test.use = "wilcox")
  markers$gene  <- rownames(markers)
  markers$clone <- cl
  markers
})
names(markers_list) <- clone_ids

all_markers <- bind_rows(markers_list)
write.csv(all_markers, file.path(OUT_DIR, "clone_markers_all.csv"), row.names = FALSE)

resistance_clone <- names(which.max(clone_prop_wide["M27", clone_ids]))
cat("\nResistance clone (highest proportion at M27):", resistance_clone, "\n")

res_markers <- markers_list[[resistance_clone]] %>%
  mutate(
    sig = case_when(
      p_val_adj < 0.05 & avg_log2FC >  0.5 ~ "Up in resistance",
      p_val_adj < 0.05 & avg_log2FC < -0.5 ~ "Down in resistance",
      TRUE ~ "NS"
    ),
    label = ifelse(sig != "NS" & rank(-abs(avg_log2FC)) <= 15, gene, "")
  )

p_volc <- ggplot(res_markers, aes(x = avg_log2FC, y = -log10(p_val_adj),
                                   colour = sig, label = label)) +
  geom_point(size = 1.2, alpha = 0.7) +
  geom_text_repel(size = 4.2, max.overlaps = 40) +
  scale_colour_manual(values = c("Up in resistance"  = "#e74c3c",
                                  "Down in resistance" = "#3498db",
                                  "NS" = "grey70")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = paste(resistance_clone, "vs other clones"),
       x = "log2 Fold Change", y = "-log10(adj. p-value)", colour = NULL) +
  theme_classic(base_size = 14) +
  theme(legend.position = "top")
ggsave(file.path(OUT_DIR, "Fig4e_volcano_resistance_clone.svg"), p_volc, width = 7, height = 5)

message("Section 5 complete.")


# =============================================================================
# SECTION 6 — ANCESTRAL STATE RECONSTRUCTION (ER MODEL)          [Results text]
# =============================================================================

message("\n=== SECTION 6: Ancestral state reconstruction ===\n")

tip_state <- setNames(
  ifelse(tree_nj$tip.label == resistance_clone, "resistant", "sensitive"),
  tree_nj$tip.label
)
tip_state["Germline"] <- "sensitive"
tip_state <- as.factor(tip_state)
print(tip_state)

fit_er <- fitMk(tree_nj, tip_state, model = "ER")
anc    <- ancr(fit_er)

state_cols <- c("sensitive" = "#3498db", "resistant" = "#e74c3c")

svg(file.path(OUT_DIR, "S_ancestral_reconstruction_ER.svg"), width = 7, height = 6)
plot(fit_er,
     tip.color = state_cols[tip_state[tree_nj$tip.label]],
     main = "Ancestral state reconstruction (ER model)")
nodelabels(pie = anc$ace, piecol = state_cols, cex = 0.6)
tiplabels(pie = to.matrix(tip_state[tree_nj$tip.label], c("sensitive", "resistant")),
          piecol = state_cols, cex = 0.4)
legend("bottomleft", legend = c("Sensitive", "Resistant"), fill = state_cols, bty = "n", cex = 0.9)
dev.off()

cat("\nPosterior probabilities at internal nodes (ER model):\n")
print(round(anc$ace, 3))

message("Section 6 complete.")


# =============================================================================
# SECTION 7 — SUMMARY
# =============================================================================

summary_df <- data.frame(
  clone     = clone_ids,
  M0_prop   = clone_prop_wide["M0", clone_ids],
  M3_prop   = clone_prop_wide["M3", clone_ids],
  M27_prop  = clone_prop_wide["M27", clone_ids],
  phenotype = tip_state[clone_ids],
  n_markers = sapply(clone_ids, function(cl)
    sum(markers_list[[cl]]$p_val_adj < 0.05 & markers_list[[cl]]$avg_log2FC > 0))
)
write.csv(summary_df, file.path(OUT_DIR, "pipeline_summary.csv"), row.names = FALSE)
print(summary_df)

message("\n=== STAGE 5 (PHYLOGENETICS) COMPLETE ===")
message("All outputs saved to: ", OUT_DIR)
