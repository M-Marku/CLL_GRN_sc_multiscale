# =============================================================================
# functions.R — functions used in more than one pipeline script.
# Source after config.R: source("config.R"); source("functions.R")
#
# Contents:
#   1. clean_aucell_colnames()  — was duplicated identically in AUCell_visualisation.R,
#                                  RSS_calculation.R, corr_aucell_fateProb.R
#   2. extract_timepoint()      — was duplicated identically in the same three scripts
#   3. calc_rss_manual()        — NEW. Replaces SCENIC::calcRSS(), which is no longer
#                                  installable in R. Same Jensen-Shannon-divergence-based
#                                  regulon specificity score as the original SCENIC package
#                                  (Suo et al. 2018; also used in pySCENIC/AUCell workflows).
#   4. gprofiler_plot()         — was defined in AUCell_visualisation.R only, but used
#                                  conceptually in RSS_calculation.R and needed (per your
#                                  note) in corr_aucell_fateProb.R too. Centralised here.
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(gprofiler2)
library(scales)


# ── 1. Clean AUCell column names ---------------------------------------------
# pySCENIC's AUCell CSV export leaves "..." in some regulon column names
# (from R's automatic name-mangling on export/import). Strips it.
clean_aucell_colnames <- function(df) {
  colnames(df) <- gsub("\\.\\.\\.", "", colnames(df))
  df
}


# ── 2. Extract timepoint from cell barcode -----------------------------------
# Cell names follow the pattern "..._M<number>_...", e.g. "CLL_CPT1_M0_AAAC...".
# Returns a factor with levels M0 < M3 < M27 (from TIMEPOINT_LEVELS in config.R),
# named by the input cell names — safe to index by barcode downstream.
extract_timepoint <- function(cell_names, levels = TIMEPOINT_LEVELS) {
  m_values <- sub(".*_M([0-9]+)_.*", "\\1", cell_names)
  tp <- factor(paste0("M", m_values), levels = levels)
  names(tp) <- cell_names
  tp
}


# ── 3. Manual Regulon Specificity Score (replaces SCENIC::calcRSS) -----------
# SCENIC (R package) is no longer installable; this reimplements its RSS
# metric directly: RSS(regulon, group) = 1 - sqrt(JSD(P, Q)), where
#   P = the regulon's activity profile across cells, normalised to sum to 1
#   Q = the indicator profile of group membership, normalised to sum to 1
#   JSD = Jensen-Shannon divergence (base-2 log)
# Reference: Suo et al. 2018, Cell Reports (regulon specificity score), as
# used in the SCENIC/pySCENIC workflows.
#
# auc_mat:      cells x regulons matrix (i.e. AUCell_df as loaded — NOT transposed)
# cell_groups:  named vector/factor, names = cell barcodes (must match rownames
#               of auc_mat), values = group label (e.g. "M0"/"M3"/"M27", or
#               "M27"/"Other" for a 2-group comparison)
#
# Returns: regulons x groups matrix of RSS values (higher = more group-specific).

.jsd <- function(p, q) {
  m <- 0.5 * (p + q)
  kl <- function(a, b) {
    idx <- a > 0
    sum(a[idx] * log2(a[idx] / b[idx]))
  }
  0.5 * kl(p, m) + 0.5 * kl(q, m)
}

calc_rss_manual <- function(auc_mat, cell_groups) {
  cells   <- intersect(rownames(auc_mat), names(cell_groups))
  auc_mat <- auc_mat[cells, , drop = FALSE]
  cell_groups <- cell_groups[cells]

  groups   <- levels(as.factor(cell_groups))
  regulons <- colnames(auc_mat)

  rss_mat <- matrix(NA_real_, nrow = length(regulons), ncol = length(groups),
                     dimnames = list(regulons, groups))

  for (g in groups) {
    q <- as.numeric(cell_groups == g)
    q <- q / sum(q)
    for (reg in regulons) {
      p <- auc_mat[, reg]
      if (sum(p) == 0) { rss_mat[reg, g] <- NA_real_; next }
      p <- p / sum(p)
      rss_mat[reg, g] <- 1 - sqrt(.jsd(p, q))
    }
  }
  rss_mat
}
# NOTE: no local test dataset was available to validate this reimplementation
# numerically against the original SCENIC::calcRSS output. If you still have
# any RSS values computed with the old package (e.g. in an old results file
# or a previous manuscript draft), it's worth a quick sanity check that the
# ranking of top regulons matches before relying on this for the paper.


# ── 4. gprofiler-based pathway enrichment plot --------------------------------
# Unchanged logic from AUCell_visualisation.R — moved here since it's called
# from more than one script (regulon-target enrichment in the RSS stage, and
# fate-associated regulon enrichment in the correlation stage).
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
