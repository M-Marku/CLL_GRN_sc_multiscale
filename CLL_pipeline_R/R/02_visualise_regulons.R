# =============================================================================
# Heatmap of top-variance AUCell regulon activity, t-SNE embedding coloured
# by timepoint and by individual regulon activity.
# =============================================================================

source("CLL_pipeline_R/config.R")
source("CLL_pipeline_R/functions.R")

library(pheatmap)
library(ggplot2)
library(tidyverse)
library(Rtsne)
library(patchwork)
library(grid)
library(viridisLite)

# ── Load AUCell results --------------------------------------------------------
AUCell_df <- read.csv(AUCELL_CSV, header = TRUE, row.names = 1, sep = "\t")
AUCell_df <- clean_aucell_colnames(AUCell_df)

cell_names <- rownames(AUCell_df)
D_values   <- extract_timepoint(cell_names)

annotation_col <- data.frame(TimePoint = D_values)
rownames(annotation_col) <- cell_names
ann_colors <- list(TimePoint = TIMEPOINT_COLORS)

# ── Top-variance TF filtering -------------------------------------
top_n       <- 100
nonzero_mat <- AUCell_df[, colSums(AUCell_df != 0) > 0]
top_genes   <- head(order(apply(nonzero_mat, 2, var), decreasing = TRUE), top_n)
filtered_mat   <- nonzero_mat[, top_genes]
TFs_ofInterest <- colnames(filtered_mat)


# ── Gap heatmap function -------------------------------------------------------
plot_gap_heatmap <- function(mat,
                              top_rows,
                              bottom_rows,
                              annotation_col,
                              ann_colors,
                              gap_label   = "···",
                              palette_n   = 100,
                              body_top    = 0.90,
                              body_bottom = 0.08,
                              rect_x      = 0.50,
                              rect_width  = 0.68,
                              label_x     = 0.755,
                              ...) {

  stopifnot(all(top_rows    %in% rownames(mat)))
  stopifnot(all(bottom_rows %in% rownames(mat)))

  top_mat    <- mat[top_rows,    , drop = FALSE]
  bottom_mat <- mat[bottom_rows, , drop = FALSE]

  gap_mat           <- matrix(colMeans(mat), nrow = 1, ncol = ncol(mat))
  rownames(gap_mat) <- ".gap."
  colnames(gap_mat) <- colnames(mat)

  display_mat <- rbind(top_mat, gap_mat, bottom_mat)
  n_total     <- nrow(display_mat)
  n_top       <- nrow(top_mat)

  row_labels            <- rownames(display_mat)
  row_labels[n_top + 1] <- gap_label
  names(row_labels)     <- rownames(display_mat)

  pal <- viridis(palette_n, option = "D")

  p <- pheatmap(
    display_mat,
    color              = pal,
    labels_row         = row_labels,
    cluster_rows       = FALSE,
    gaps_row           = c(n_top, n_top + 1),
    annotation_col     = annotation_col,
    annotation_colors  = ann_colors,
    show_colnames      = FALSE,
    ...
  )

  body_height <- body_top - body_bottom
  row_h       <- body_height / n_total
  gap_top_y   <- body_top    - n_top * row_h
  gap_bot_y   <- gap_top_y   - row_h
  gap_mid_y   <- (gap_top_y  + gap_bot_y) / 2

  pushViewport(viewport(
    x = unit(0, "npc"), y = unit(0, "npc"),
    width = unit(1, "npc"), height = unit(1, "npc"),
    just = c("left", "bottom")
  ))

  grid.rect(
    x      = unit(rect_x,     "npc"),
    y      = unit(gap_mid_y,  "npc"),
    width  = unit(rect_width, "npc"),
    height = unit(row_h + 0.004, "npc"),
    gp     = gpar(fill = "white", col = NA)
  )

  grid.text(
    label = gap_label,
    x     = unit(label_x, "npc"),
    y     = unit(gap_mid_y, "npc"),
    gp    = gpar(col = "grey50", fontsize = 10, fontface = "bold")
  )

  popViewport()
  invisible(p)
}

# ── Define which TFs to show in the gap heatmap --------------------------------
top_show <- c("TAF7", "MSI2", "ENO1", "IRF8", "SPIB", "YBX1", "HMGB1",
              "ARID5B", "ILF2", "ZBTB7A", "KLF6", "ZNF160", "CREB3L2")
bot_show <- c("ETV5", "JUN", "JUND", "DDIT3", "ATF4", "JUNB", "TCF4",
              "MEF2C", "BCL11A", "ATF6B", "FOXP1", "IRF1", "STAT1")

svg(file.path(PLOTS_DIR, "heatmap_gap_TFs_AUCell.svg"), width = 9, height = 6)
plot_gap_heatmap(
  mat            = t(filtered_mat),
  top_rows       = top_show,
  bottom_rows    = bot_show,
  annotation_col = annotation_col,
  ann_colors     = ann_colors,
  gap_label      = "···",
  fontsize       = 10,
  main           = paste("Top", top_n, "TFs regulon activity")
)
dev.off()


# ── t-SNE on regulon activity ---------------------------------------------------
aucMat     <- AUCell_df[, apply(AUCell_df, 2, sd) > 0]
auc_scaled <- scale(aucMat)
pca        <- prcomp(auc_scaled, center = FALSE, scale. = FALSE)
n_pcs      <- 30

set.seed(100)
tsne <- Rtsne(pca$x[, 1:n_pcs], perplexity = 30, theta = 0.5, verbose = TRUE)

tsne_df <- data.frame(
  tSNE1 = tsne$Y[, 1],
  tSNE2 = tsne$Y[, 2],
  Time  = D_values,
  Cell  = rownames(AUCell_df)
)

tsne_time <- ggplot(tsne_df, aes(tSNE1, tSNE2, color = Time)) +
  geom_point(size = 0.8, alpha = 0.7) +
  theme_classic() +
  theme(legend.title = element_text(size = 14),
        legend.text  = element_text(size = 12)) +
  labs(title = "Regulon-based t-SNE", color = "Time point") +
  scale_color_manual(values = TIMEPOINT_COLORS)

ggsave(file.path(PLOTS_DIR, "tSNE_regulons_time_annot.svg"),
       plot = tsne_time, width = 5, height = 4, dpi = 300)

regulon_names <- c("ARID5B", "ATF4", "DDIT3", "ETV5", "FOS",
                    "HMGB1", "TAF7", "JUN", "JUNB", "JUND")

plot_tf_tsne_grid <- function(tsne_df, aucMat, tfs,
                               ncol = 5, point_size = 0.7, show_legend = TRUE) {
  plots <- list()
  for (tf in tfs) {
    if (!tf %in% colnames(aucMat)) { warning(paste(tf, "not found, skipping")); next }
    tsne_df$activity <- aucMat[tsne_df$Cell, tf]
    p <- ggplot(tsne_df, aes(tSNE1, tSNE2, color = activity)) +
      geom_point(size = point_size) +
      scale_color_viridis_c(option = "D") +
      theme_classic() +
      labs(title = tf, color = "Regulon activity") +
      theme(plot.title      = element_text(hjust = 0.5, size = 14),
            legend.text     = element_text(size = 12),
            legend.position = ifelse(show_legend, "right", "none"))
    plots[[tf]] <- p
  }
  wrap_plots(plots, ncol = ncol)
}

plots_arranged <- plot_tf_tsne_grid(tsne_df, aucMat, regulon_names)
ggsave(file.path(PLOTS_DIR, "tSNE_regulons_top10M27.svg"),
       plot = plots_arranged, width = 24, height = 10, dpi = 300)

message("Stage 2 complete: AUCell heatmap and t-SNE plots saved to ", PLOTS_DIR)

saveRDS(TFs_ofInterest, file.path(RESULTS_DIR, "top100_variance_TFs.rds"))
