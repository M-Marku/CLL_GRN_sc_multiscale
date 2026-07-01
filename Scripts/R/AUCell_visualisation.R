library(pheatmap)
library(ggplot2)
library(tidyverse)
library(gprofiler2)
library(scales)
library(ggpubr)
library(Rtsne)
library(patchwork)
library(grid)
library(viridisLite)

# Load the AUCell results 
AUCell_df <- read.csv("Results/CLL_only/pyScenic_results/AUcell_CLL_only.csv", 
                      header = T, row.names = 1, sep = "\t")

colnames(AUCell_df) <- gsub("\\.\\.\\.", "", colnames(AUCell_df))

cell_names <- rownames(AUCell_df)
M_values   <- sub(".*_M([0-9]+)_.*", "\\1", cell_names)
D_values   <- factor(paste0("M", M_values), levels = c("M0", "M3", "M27"))

annotation_col <- data.frame(TimePoint = D_values)
rownames(annotation_col) <- cell_names

category_colors <- c("M0" = "#66C2A5", "M3" = "#FC8D62", "M27" = "#8DA0CB")
ann_colors      <- list(TimePoint = category_colors)

# ── Top TF filtering (unchanged) ──────────────────────────────────────────────
top_n       <- 100
nonzero_mat <- AUCell_df[, colSums(AUCell_df != 0) > 0]
top_genes   <- head(order(apply(nonzero_mat, 2, var), decreasing = TRUE), top_n)
filtered_mat <- nonzero_mat[, top_genes]
TFs_ofInterest <- colnames(filtered_mat)

# ── Gap heatmap function ───────────────────────────────────────────────────────
# Rows of the heatmap = TFs (we transpose filtered_mat).
# top_rows / bottom_rows are TF names (character vectors).
# Everything between them is replaced by a single "···" spacer row.

plot_gap_heatmap <- function(mat,
                             top_rows,
                             bottom_rows,
                             annotation_col,
                             ann_colors,
                             gap_label     = "···",
                             palette_n     = 100,
                             body_top      = 0.90,
                             body_bottom   = 0.08,
                             rect_x        = 0.50,
                             rect_width    = 0.68,
                             label_x       = 0.755,
                             ...) {
  
  stopifnot(all(top_rows    %in% rownames(mat)))
  stopifnot(all(bottom_rows %in% rownames(mat)))
  
  # ── 1. Build display matrix ──────────────────────────────────────────────────
  top_mat    <- mat[top_rows,    , drop = FALSE]
  bottom_mat <- mat[bottom_rows, , drop = FALSE]
  
  # Gap row: column means → renders as mid-colour (will be whited-out anyway)
  gap_mat           <- matrix(colMeans(mat), nrow = 1, ncol = ncol(mat))
  rownames(gap_mat) <- ".gap."
  colnames(gap_mat) <- colnames(mat)
  
  display_mat <- rbind(top_mat, gap_mat, bottom_mat)
  n_total     <- nrow(display_mat)
  n_top       <- nrow(top_mat)
  
  # ── 2. Row labels: real names + gap_label for the spacer ─────────────────────
  row_labels        <- rownames(display_mat)
  row_labels[n_top + 1] <- gap_label
  names(row_labels) <- rownames(display_mat)
  
  # ── 3. Viridis palette ───────────────────────────────────────────────────────
  pal <- viridis(palette_n, option = "D")
  
  # ── 4. Draw pheatmap ─────────────────────────────────────────────────────────
  p <- pheatmap(
    display_mat,
    color          = pal,
    labels_row     = row_labels,
    cluster_rows   = FALSE,
    gaps_row       = c(n_top, n_top + 1),   # thin lines above and below gap row
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    show_colnames  = FALSE,
    ...
  )
  
  # ── 5. Grid overlay: white-out the gap row, redraw "···" ─────────────────────
  body_height <- body_top - body_bottom
  row_h       <- body_height / n_total
  gap_top_y   <- body_top    - n_top       * row_h
  gap_bot_y   <- gap_top_y   - row_h
  gap_mid_y   <- (gap_top_y  + gap_bot_y)  / 2
  
  pushViewport(viewport(
    x = unit(0, "npc"), y = unit(0, "npc"),
    width = unit(1, "npc"), height = unit(1, "npc"),
    just = c("left", "bottom")
  ))
  
  # White rectangle over the gap row cells
  grid.rect(
    x      = unit(rect_x,     "npc"),
    y      = unit(gap_mid_y,  "npc"),
    width  = unit(rect_width, "npc"),
    height = unit(row_h + 0.004, "npc"),
    gp     = gpar(fill = "white", col = NA)
  )
  
  # "···" label in the row-label margin
  grid.text(
    label = gap_label,
    x     = unit(label_x, "npc"),
    y     = unit(gap_mid_y, "npc"),
    gp    = gpar(col = "grey50", fontsize = 10, fontface = "bold")
  )
  
  popViewport()
  invisible(p)
}

# ── Define which TFs to show ───────────────────────────────────────────────────
# Option A: top N by variance (already ranked in filtered_mat columns)
#   Show first 15, skip the middle, show last 10.
all_tfs  <- colnames(filtered_mat)   # already ranked by variance
top_tfs <- c("TAF7", "MSI2", "ENO1", "IRF8", "SPIB", "YBX1", "HMGB1", "ARID5B", "ILF2", "ZBTB7A", "KLF6", "ZNF160", "CREB3L2")
bottom_tfs <- c("ETV5", "JUN", "JUND", "DDIT3", "ATF4", "JUNB", "TCF4", "MEF2C", "BCL11A", "ATF6B", "FOXP1", "IRF1", "STAT1")
#top_show <- all_tfs[1:15]
#bot_show <- tail(all_tfs, 10)

# Option B: hand-pick biologically meaningful groups, e.g.:
# top_show <- c("JUN", "JUND", "JUNB", "ATF4", "DDIT3", "ETV5", "KLF6",
#               "TAF7", "ARID5B", "HMGB1")
# bot_show <- c("POU2F2", "SPIB", "IRF4", "BATF", "NFKB1", "RELA",
#               "IRF3", "MYC", "STAT1", "NRF1")

# ── Save the gap heatmap ───────────────────────────────────────────────────────
svg("Plots/CLL_only/heatmap_gap_TFs_AUCell.svg", width = 9, height = 6)

plot_gap_heatmap(
  mat            = t(filtered_mat),   # rows = TFs, cols = cells
  top_rows       = top_show,
  bottom_rows    = bot_show,
  annotation_col = annotation_col,
  ann_colors     = ann_colors,
  gap_label      = "···",
  # pheatmap pass-through args:
  #cutree_cols    = 3,
  fontsize       = 10,
  main           = paste("Top", top_n, "TFs regulon activity")
)

dev.off()

# ── t-SNE (unchanged below) ───────────────────────────────────────────────────
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
  scale_color_manual(values = c("M0" = "#66C2A5", "M3" = "#FC8D62", "M27" = "#8DA0CB"))

ggsave("Plots/CLL_only/tSNE_regulons_time_annot.svg",
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
      theme(plot.title    = element_text(hjust = 0.5, size = 14),
            legend.text   = element_text(size = 12),
            legend.position = ifelse(show_legend, "right", "none"))
    plots[[tf]] <- p
  }
  wrap_plots(plots, ncol = ncol)
}

plots_arranged <- plot_tf_tsne_grid(tsne_df, aucMat, regulon_names)
ggsave("Plots/CLL_only/tSNE_regulons_top10M27.svg",
       plot = plots_arranged, width = 24, height = 10, dpi = 300)

# ── Pathway enrichment (unchanged) ───────────────────────────────────────────
gprofiler_plot <- function(gene_list, custom.bg, gprof_sources,
                           max_set_size, min_set_size, correction_method,
                           n_terms_per_db, title) {
  gostres <- gost(gene_list, organism = "hsapiens", ordered_query = FALSE,
                  evcodes = TRUE, correction_method = correction_method,
                  domain_scope = "annotated", custom_bg = custom.bg,
                  sources = gprof_sources, exclude_iea = FALSE)
  if (is.null(gostres)) {
    message("No significant enrichment results to display")
    return(invisible(NULL))
  }
  gostres.df <- gostres$result
  gostres.df <- dplyr::filter(gostres.df,
                              term_size <= max_set_size & term_size > min_set_size)
  dup <- duplicated(gostres.df$term_name)
  gostres.df$term_name[dup] <- paste(gostres.df$source[dup],
                                     gostres.df$term_name[dup], sep = "_")
  tmp <- gostres.df %>%
    arrange(p_value) %>%
    group_by(source) %>%
    slice(1:n_terms_per_db) %>%
    as.data.frame()
  tmp$term_name <- factor(tmp$term_name,
                          levels = fct_reorder(rev(tmp$term_name),
                                               rev(-log10(tmp$p_value))))
  ggplot(tmp, aes(-log10(p_value), term_name, color = source, size = recall)) +
    geom_point() +
    theme_bw() +
    theme(legend.position = "bottom",
          axis.title.y    = element_blank(),
          axis.text.y     = element_text(size = 14)) +
    labs(x = "-log10(p-value)", size = "Ratio", title = title) +
    guides(colour = guide_legend(nrow = 1, order = 1, direction = "vertical",
                                 title.position = "top"),
           size   = guide_legend(nrow = 1, order = 2, direction = "vertical",
                                 title.position = "top")) +
    scale_y_discrete(labels = label_wrap(width = 30))
}

Guniverse <- decoupleR::get_collectri(organism = "human", split_complexes = FALSE)
TF_universe <- Guniverse[Guniverse$source %in% TFs_ofInterest, ]
top_TFs     <- unique(c(TF_universe$source, TF_universe$target))

pathways_plot <- gprofiler_plot(
  top_TFs, custom.bg = NULL,
  gprof_sources = c("REAC", "KEGG", "GO:BP"),
  max_set_size = 500, min_set_size = 10,
  correction_method = "fdr", n_terms_per_db = 10,
  title = "Pathway enrichment of top 10 M27 TFs and their targets"
)
ggsave("Plots/CLL_only/pathway_enrichment_top_M3_driver_genes.svg",
       plot = pathways_plot, width = 6.5, height = 6, dpi = 300)
