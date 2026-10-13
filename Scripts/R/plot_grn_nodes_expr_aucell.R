# Two-row figure: expression (top) and AUCell regulon score (bottom) of GRN / Boolean-model
# nodes across timepoints. AUCell exists only for TFs with a pySCENIC regulon; nodes without
# one (kinases, BCL2 family...) get an empty panel in the bottom row so columns stay aligned.

library(Seurat)
library(tidyverse)
library(patchwork)

# ---- edit ---------------------------------------------------------------------------------
SEU_PATH    <- "Data/CLL_only/cll_obj_tfAct_AUCell.rds"
TIME_COL    <- "time_point"
AUCELL_ASSAY <- "AUCell"
OUT_FILE    <- "Plots/CLL_only/GRN_nodes_expr_aucell.svg"
NODES       <- c("BTK", "PLCG2", "MTOR", "AKT1", "NFKB1", "BCL2L1", "BCL2", "MCL1", "JUN", "MYC")
TP_LEVELS   <- c("M0", "M3", "M27")
TP_COLORS   <- c(M0 = "#66C2A5", M3 = "#FC8D62", M27 = "#8DA0CB")
# --------------------------------------------------------------------------------------------

seu <- readRDS(SEU_PATH)
tp  <- setNames(factor(seu[[TIME_COL]][, 1], levels = TP_LEVELS), colnames(seu))

# long data frame (cell x node) from one assay
get_long <- function(assay) {
  mat   <- GetAssayData(seu, assay = assay, layer = "data")
  found <- intersect(NODES, rownames(mat))
  if (length(setdiff(NODES, found)))
    message(assay, ": not found -> ", paste(setdiff(NODES, found), collapse = ", "))
  as.data.frame(t(as.matrix(mat[found, , drop = FALSE]))) %>%
    mutate(time = tp[rownames(.)]) %>%
    pivot_longer(-time, names_to = "node", values_to = "value") %>%
    mutate(node = factor(node, levels = NODES))     # keep all nodes so panels align
}

plot_row <- function(df, ylab) {
  ggplot(df, aes(time, value, fill = time)) +
    geom_violin(scale = "width", linewidth = 0.3, alpha = 0.85) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", linewidth = 0.3) +
    facet_wrap(~ node, nrow = 1, scales = "free_y", drop = FALSE) +
    scale_fill_manual(values = TP_COLORS) +
    labs(x = NULL, y = ylab) +
    theme_classic(base_size = 10) +
    theme(legend.position = "none", strip.background = element_blank(),
          strip.text = element_text(face = "bold"))
}

p <- plot_row(get_long("RNA"), "Expression (log-normalised)") /
     plot_row(get_long(AUCELL_ASSAY), "AUCell regulon score")

dir.create(dirname(OUT_FILE), showWarnings = FALSE, recursive = TRUE)
ggsave(OUT_FILE, p, width = max(6, 1.5 * length(NODES)), height = 5)
