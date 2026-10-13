# =============================================================================
# Spearman correlation between regulon (AUCell) activity and CellRank fate
# probabilities, per lineage (M3, M27); heatmap of top correlated regulons;
# pathway enrichment of the fate-associated regulons and their targets.
# =============================================================================

source("Scripts/config.R")
source("Scripts/functions.R")

library(tidyverse)
library(pheatmap)
library(viridis)

# ── Load AUCell + fate probabilities -------------------------------------------
AUCell_df <- read.csv(AUCELL_CSV, header = TRUE, row.names = 1, sep = "\t")
AUCell_df <- clean_aucell_colnames(AUCell_df)

fate_probs <- read.csv(FATE_PROBS_CSV, header = TRUE, row.names = 1)
colnames(fate_probs) <- c("M3 lineage", "M27 lineage")

common_cells <- intersect(rownames(AUCell_df), rownames(fate_probs))
fate_probs   <- fate_probs[common_cells, ]
AUCell_df    <- AUCell_df[common_cells, ]

# ── Spearman correlation: regulon activity vs. fate probability ----------------
cor_matrix <- matrix(
  NA, nrow = ncol(AUCell_df), ncol = ncol(fate_probs),
  dimnames = list(colnames(AUCell_df), colnames(fate_probs))
)
pval_matrix <- cor_matrix

for (fate in colnames(fate_probs)) {
  for (reg in colnames(AUCell_df)) {
    res <- cor.test(fate_probs[[fate]], AUCell_df[[reg]], method = "spearman")
    cor_matrix[reg, fate]  <- res$estimate
    pval_matrix[reg, fate] <- res$p.value
  }
}
pval_adj <- apply(pval_matrix, 2, p.adjust, method = "BH")

fateM3_regulons <- tibble(
  regulon     = rownames(cor_matrix),
  correlation = cor_matrix[, "M3 lineage"],
  p_adj       = pval_adj[, "M3 lineage"]
) %>% filter(p_adj < 0.05) %>% arrange(desc(abs(correlation)))

fateM27_regulons <- tibble(
  regulon     = rownames(cor_matrix),
  correlation = cor_matrix[, "M27 lineage"],
  p_adj       = pval_adj[, "M27 lineage"]
) %>% filter(p_adj < 0.05) %>% arrange(desc(abs(correlation)))

write.csv(fateM3_regulons,  file.path(RESULTS_DIR, "fate_regulons_M3.csv"),  row.names = FALSE)
write.csv(fateM27_regulons, file.path(RESULTS_DIR, "fate_regulons_M27.csv"), row.names = FALSE)

# ── Heatmap of top regulons per lineage (unchanged figure) ---------------------
top_regs <- unique(c(
  fateM3_regulons$regulon[1:20],
  fateM27_regulons$regulon[1:20]
))
heatmap_mat <- cor_matrix[top_regs, ]

svg(file.path(PLOTS_DIR, "corr_AUCell_FateProb_t.svg"), width = 8, height = 3)
pheatmap(
  t(heatmap_mat),
  color = viridis(100),
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  border_color = NA,
  main = "Correlation between regulon activity and fate probabilities"
)
dev.off()

message("Correlation analysis complete.")


# =============================================================================
# SECTION B — Pathway enrichment of fate-associated regulons (per lineage)
# Uses ALL significant regulons per lineage (p_adj < 0.05), not an arbitrary
# top-N — the significance filter above already does that job.
# =============================================================================

enrich_fate_lineage <- function(regulon_df, lineage_label, out_filename) {
  if (nrow(regulon_df) == 0) {
    message("No significant regulons for ", lineage_label, " — skipping enrichment.")
    return(invisible(NULL))
  }
  Guniverse   <- decoupleR::get_collectri(organism = "human", split_complexes = FALSE)
  TF_universe <- Guniverse[Guniverse$source %in% regulon_df$regulon, ]
  targets     <- unique(c(TF_universe$source, TF_universe$target))

  p <- gprofiler_plot(
    targets, custom.bg = NULL,
    gprof_sources = c("REAC", "KEGG", "GO:BP"),
    max_set_size = 500, min_set_size = 10,
    correction_method = "fdr", n_terms_per_db = 10,
    title = paste("Pathway enrichment —", lineage_label, "fate-associated regulons")
  )
  if (!is.null(p)) {
    # NOTE: no prior figure/size existed for this analysis (it's newly added) —
    # picked a size consistent with the other enrichment plot in this pipeline
    # (pathway_enrichment_top_M3_driver_genes.svg, 6.5 x 6); adjust freely.
    ggsave(file.path(PLOTS_DIR, out_filename), plot = p, width = 6.5, height = 6, dpi = 300)
  }
  invisible(p)
}

enrich_fate_lineage(fateM3_regulons,  "M3",  "pathway_enrichment_M3_fate_regulons.svg")
enrich_fate_lineage(fateM27_regulons, "M27", "pathway_enrichment_M27_fate_regulons.svg")

message("Stage 4 complete: correlation heatmap + fate-lineage pathway enrichment saved.")
