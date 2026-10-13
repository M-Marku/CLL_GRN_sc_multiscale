# =============================================================================
# Regulon Specificity Score (RSS) for M27 vs. rest, ranking plot, violin plots
# of top regulons, Wilcoxon test for M27-up/down regulons, and pathway
# enrichment of the M27-associated regulons.
# =============================================================================

source("Scripts/config.R")
source("Scripts/functions.R")

library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)

# ── Load AUCell + assign timepoints -------------------------------------------
AUCell_df <- read.csv(AUCELL_CSV, header = TRUE, row.names = 1, sep = "\t")
AUCell_df <- clean_aucell_colnames(AUCell_df)
D_values  <- extract_timepoint(rownames(AUCell_df))

# ── Wilcoxon test: M27 vs. M0+M3 for every regulon -----------------------------
df_long <- AUCell_df %>%
  as.data.frame() %>%
  mutate(Cell = rownames(.), Group = D_values[Cell]) %>%
  pivot_longer(cols = -c(Cell, Group), names_to = "Regulon", values_to = "Activity")

wilcox_results <- df_long %>%
  group_by(Regulon) %>%
  summarise(
    p_value = wilcox.test(
      Activity[Group == "M27"],
      Activity[Group %in% c("M0", "M3")]
    )$p.value,
    direction = median(Activity[Group == "M27"]) -
      median(Activity[Group %in% c("M0", "M3")]),
    .groups = "drop"
  ) %>%
  mutate(FDR = p.adjust(p_value, method = "fdr")) %>%
  arrange(FDR)

sig_up_M27   <- wilcox_results %>% filter(FDR < 0.05 & direction > 0 & p_value < 0.05)
sig_down_M27 <- wilcox_results %>% filter(FDR < 0.05 & direction < 0 & p_value < 0.05)

write.csv(wilcox_results, file.path(RESULTS_DIR, "wilcoxon_M27_vs_rest.csv"), row.names = FALSE)
cat("Significant M27-up regulons:", nrow(sig_up_M27), "\n")
cat("Significant M27-down regulons:", nrow(sig_down_M27), "\n")

# ── RSS: M27 vs. Other (2-group), using the manual reimplementation ------------
group_2 <- factor(ifelse(D_values == "M27", "M27", "Other"), levels = c("M27", "Other"))
names(group_2) <- names(D_values)

rss_2groups <- calc_rss_manual(as.matrix(AUCell_df), group_2)

rss_df <- data.frame(
  regulon = rownames(rss_2groups),
  RSS     = rss_2groups[, "M27"]
) %>%
  arrange(desc(RSS)) %>%
  mutate(rank = row_number())

write.csv(rss_df, file.path(RESULTS_DIR, "RSS_M27_vs_rest.csv"), row.names = FALSE)

# Label the TFs found significantly up-regulated in M27 by the Wilcoxon test
tfs_to_label <- sig_up_M27$Regulon

RSS_plot <- ggplot(rss_df, aes(x = rank, y = RSS)) +
  geom_point(color = "steelblue", size = 2) +
  geom_point(data = subset(rss_df, regulon %in% tfs_to_label),
             color = "red", size = 3) +
  geom_text_repel(
    data = subset(rss_df, regulon %in% tfs_to_label),
    aes(label = regulon),
    color = "red", size = 4, nudge_x = 0.5, nudge_y = 0.01,
    max.overlaps = Inf, box.padding = 0.4, force = 2,
    segment.color = "grey50", direction = "y", seed = 42
  ) +
  labs(x = "Regulon", y = "RSS", title = "RSS score, M27") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(), axis.text.x = element_blank())

ggsave(file.path(PLOTS_DIR, "RSS_plot_M27.svg"), plot = RSS_plot,
       device = "svg", width = 5, height = 7)

# ── Violin plots of top regulons by RSS ----------------------------------------
n_TFs_of_interest <- 9
top_TFs <- rss_df %>% slice_head(n = n_TFs_of_interest) %>% pull(regulon)

df_plot <- AUCell_df[, top_TFs, drop = FALSE] %>%
  as.data.frame() %>%
  mutate(Cell = rownames(.), Group = D_values[Cell]) %>%
  pivot_longer(cols = all_of(top_TFs), names_to = "Regulon", values_to = "Activity")

plot_regulon_distribution <- ggplot(df_plot, aes(x = Group, y = Activity, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.5) +
  facet_wrap(~ Regulon, scales = "free_y", nrow = 4) +
  theme_bw(base_size = 14) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Cell category", y = "Regulon activity",
       title = paste0("Top ", n_TFs_of_interest, " regulons in M27 category")) +
  theme(panel.grid = element_blank(), strip.text = element_text(size = 12, face = "bold"))

ggsave(file.path(PLOTS_DIR, "top_9_regulon_AUCell.svg"),
       plot = plot_regulon_distribution, width = 8, height = 6, dpi = 300)

message("RSS + Wilcoxon analysis complete.")


# =============================================================================
# SECTION A — ORIGINAL pathway enrichment: top-10 M27 TFs + CollecTRI targets,
#             gprofiler ORA.
# =============================================================================

TFs_ofInterest <- readRDS(file.path(RESULTS_DIR, "top100_variance_TFs.rds"))

Guniverse    <- decoupleR::get_collectri(organism = "human", split_complexes = FALSE)
TF_universe  <- Guniverse[Guniverse$source %in% TFs_ofInterest, ]
top_TFs_targets <- unique(c(TF_universe$source, TF_universe$target))

pathways_plot_topN <- gprofiler_plot(
  top_TFs_targets, custom.bg = NULL,
  gprof_sources = c("REAC", "KEGG", "GO:BP"),
  max_set_size = 500, min_set_size = 10,
  correction_method = "fdr", n_terms_per_db = 10,
  title = "Pathway enrichment of top 10 M27 TFs and their targets"
)
ggsave(file.path(PLOTS_DIR, "pathway_enrichment_top_M3_driver_genes.svg"),
       plot = pathways_plot_topN, width = 6.5, height = 6, dpi = 300)


# =============================================================================
# SECTION B — Ranked-based (preranked GSEA) alternative, no fixed cutoff.
#             Ranks ALL regulons by RSS and tests whether pathway gene sets are
#             concentrated toward the M27-specific end of that ranking.
# =============================================================================

library(fgsea)
library(msigdbr)

ranked_tfs <- setNames(rss_df$RSS, rss_df$regulon)
ranked_tfs <- sort(ranked_tfs, decreasing = TRUE)

# Combine Reactome + KEGG + GO:BP to mirror the original gprof_sources choice
msig_sets <- bind_rows(
  msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME"),
  msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
)
pathway_list <- split(msig_sets$gene_symbol, msig_sets$gs_name)

set.seed(42)
gsea_res <- fgsea(
  pathways = pathway_list,
  stats    = ranked_tfs,
  minSize  = 3,     # small on purpose — ranking universe is ~568 TFs, not the whole transcriptome
  maxSize  = 500,
  eps      = 0
)
gsea_res <- gsea_res[order(gsea_res$padj), ]

write.csv(as.data.frame(gsea_res[, c("pathway", "NES", "pval", "padj", "size")]),
          file.path(RESULTS_DIR, "GSEA_ranked_M27_regulons.csv"), row.names = FALSE)

message("Stage 3 complete: RSS, Wilcoxon test, and both enrichment analyses (A: top-10 ORA, B: ranked GSEA) done.")
