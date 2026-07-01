# Load the libraries
library(tidyverse)
library(pheatmap)
library(viridis)

# Load the AUCell matrix
AUCell_df <- read.csv("Results/CLL_only/AUcell_CLL_only.csv", header = T, row.names = 1, sep = "\t")

# Remove the ... in the column names
colnames(AUCell_df) <- gsub("\\.\\.\\.", "", colnames(AUCell_df))
cell_names <- rownames(AUCell_df)

# Load the fate probabilities matrix
fate_probs <- read.csv("Results/CLL_only/CellRank_results/cell_matrix_fate_probabilities_light.csv", header = T, row.names = 1)
colnames(fate_probs) <- c("M3 lineage", "M27 lineage")

common_cells <- intersect(rownames(AUCell_df), rownames(fate_probs))

fate_probs  <- fate_probs[common_cells, ]
AUCell_df <- AUCell_df[common_cells, ]

# Estimate the Spearman correlation between regulon activity and fate probability
cor_matrix <- matrix(
  NA,
  nrow = ncol(AUCell_df),
  ncol = ncol(fate_probs),
  dimnames = list(colnames(AUCell_df), colnames(fate_probs))
)

pval_matrix <- cor_matrix

for (fate in colnames(fate_probs)) {
  for (reg in colnames(AUCell_df)) {
    res <- cor.test(
      fate_probs[[fate]],
      AUCell_df[[reg]],
      method = "spearman"
    )
    cor_matrix[reg, fate]  <- res$estimate
    pval_matrix[reg, fate] <- res$p.value
  }
}
pval_adj <- apply(pval_matrix, 2, p.adjust, method = "BH")

# Select the regulons for each lineage
fateM3_regulons <- tibble(
  regulon = rownames(cor_matrix),
  correlation = cor_matrix[, "M3 lineage"],
  p_adj = pval_adj[, "M3 lineage"]
) %>%
  filter(p_adj < 0.05) %>%
  arrange(desc(abs(correlation)))

fateM27_regulons <- tibble(
  regulon = rownames(cor_matrix),
  correlation = cor_matrix[, "M27 lineage"],
  p_adj = pval_adj[, "M27 lineage"]
) %>%
  filter(p_adj < 0.05) %>%
  arrange(desc(abs(correlation)))

# Select the top regulons per lineage and visualise 
top_regs <- unique(c(
  fateM3_regulons$regulon[1:20],
  fateM27_regulons$regulon[1:20]
))

heatmap_mat <- cor_matrix[top_regs, ]
svg("Plots/CLL_only/corr_AUCell_FateProb_t.svg", width = 8, height = 3)
pheatmap(
  t(heatmap_mat),
  color = viridis(100),
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  border_color = NA,
  main = "Correlation between regulon activity and fate probabilities"
)
dev.off()
