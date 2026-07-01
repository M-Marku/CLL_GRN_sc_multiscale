
library(Matrix)
library(RandomWalkRestartMH)
library(igraph)
library(ggraph)
library(precrec)
library(tidyverse)

# Example taken from here: https://alberto-valdeolivas.github.io/RandomWalkRestartMH/articles/RandomWalkRestartMH.html

M0_GRN <- read.csv("Results/CLL_only/pyScenic_results/SCENIC_GRN_edges_M0.csv", header = T)
network_graph <- graph_from_data_frame(M0_GRN, directed = TRUE)

# Extract only the giant component 
components <- components(network_graph, mode = "weak")
giant_component_ids <- which(components$membership == which.max(components$csize))
giant_g <- induced_subgraph(network_graph, giant_component_ids)

M0_MultiplexObject <- create.multiplex(list(GRN=giant_g))
AdjMatrix_M0 <- compute.adjacency.matrix(M0_MultiplexObject)
AdjMatrixNorm_M0 <- normalize.multiplex.adjacency(AdjMatrix_M0)
AdjMatrixNorm_M0[is.na(AdjMatrixNorm_M0)] <- 0
AdjMatrixNorm_M0[is.nan(AdjMatrixNorm_M0)] <- 0

# SeedGene <- c("BTK", 
#               "TEC", 
#               "ITK", 
#               "BMX", 
#               "TXK", 
#               "EGFR", 
#               "ERBB2",
#               "ERBB4",
#               "BLK", 
#               "JAK3",
#               "CSK",
#               "BRK")

ibrutinib_targets <- read.csv("Data/CLL_only/drug_interaction_ibrutinib.tsv", header = T, sep = "\t")
SeedGene <- ibrutinib_targets$gene
#seeds_to_test <- c("BTK", "TEC", "RELA", "NFKB1", "JUN")
#GRN_nodes <- unique(c(M0_GRN$source, M0_GRN$target))
GRN_nodes <- V(giant_g)$name
SeedGene <- intersect(SeedGene, GRN_nodes)
## We launch the algorithm with the default parameters (See details on manual)
RWR_M0_Results <- Random.Walk.Restart.Multiplex(AdjMatrixNorm_M0,
                                                 M0_MultiplexObject,
                                                SeedGene,
                                                 r = 0.7) # the higher the r, the more the walker stays around the seed

# TopResults_M0 <-
#   create.multiplexNetwork.topResults(RWR_M0_Results,
#                                      M0_MultiplexObject,
#                                      k=15)
# ## Plot the interactions of the seed
# par(mar=c(0.1,0.1,0.1,0.1))
# plot(TopResults_M0, 
#      vertex.label.color="black",
#      vertex.frame.color="#ffffff",
#      vertex.size= 20, 
#      edge.curved=.2,
#      vertex.color = ifelse(igraph::V(TopResults_M0)$name == SeedGene,"yellow",
#                            "#00CCFF"), edge.color="blue",edge.width=0.8)

## Visualise the graph with the gene scores 
prop_scores <- as.data.frame(RWR_M0_Results$RWRM_Results)
prop_scores$Score <- prop_scores$Score / sum(prop_scores$Score)
top_50_genes <- prop_scores %>%
  arrange(desc(Score)) %>%  # Sort by score descending
  slice_head(n = 200) 
top_50_genes <- top_50_genes$NodeNames

# Induce a subnetwork of these genes
sub_g <- induced_subgraph(network_graph, vids = top_50_genes)

# Plot with ggraph
top10_genes <- prop_scores$NodeNames[
  order(prop_scores$Score, decreasing = TRUE)[1:10]]

RWR_graph <- ggraph(sub_g, layout = "nicely") +
  geom_edge_link(aes(alpha = weight),
                 arrow = arrow(length = unit(1, 'mm'))) +
  geom_node_point(aes(
    size = prop_scores$Score[match(name, prop_scores$NodeNames)],
    color = prop_scores$Score[match(name, prop_scores$NodeNames)],
    stroke = ifelse(name %in% SeedGene, 2, 0.5)
  ), shape = 21, fill = "white") + 
  geom_node_text(
    aes(label = name),
    color = ifelse(V(sub_g)$name %in% top10_genes, "blue", "black"),
    repel = TRUE
  ) +
  scale_color_distiller(palette = "YlOrRd", direction = 1) +
  theme_graph() + 
  labs(
    title = "Ibrutinib Network Propagation",
    color = "Regulatory Influence",
    size = "Regulatory Influence"
  )
ggsave(RWR_graph, filename = "Plots/CLL_only/RWR_M27.svg", device = "svg", width = 15, height = 10, dpi = 300)

# Define the M27-specific TFs
#M27_TFs <- c("ARID5B", "ATF4", "DDIT3", "ETV5", "FOS", "HMGB1", "TAF7", "JUN", "JUNB", "JUND")
M27_genes <- read.table("Results/CLL_only/genes_M27_cond.txt")
M3_genes <- read.table("Results/CLL_only/genes_M3_cond.txt")

# Estimate the AUROC and AUPRC in the complete network

prop_scores$is_m27 <- ifelse(prop_scores$NodeNames %in% M27_genes$V1, 1, 0)
prop_scores$cum_m27 <- cumsum(prop_scores$is_m27)
prop_scores <- prop_scores[order(prop_scores$Score, decreasing = TRUE), ]
prop_scores$rank <- 1:nrow(prop_scores)

ggplot(prop_scores, aes(x = rank, y = cum_m27)) +
  geom_line(color = "firebrick", size = 1) +
  geom_abline(intercept = 0, slope = sum(prop_scores$is_m27)/nrow(prop_scores), linetype = "dashed") +
  labs(title = "M27 relapse enrichment curve",
       subtitle = "Does the drug perturbation (M0) point toward the relapse (M27) TFs?",
       x = "Ranked genes (by RWR score)",
       y = "Cumulative M27 gemes found") +
  theme_minimal()

precrec_obj <- evalmod(scores = prop_scores$Score, labels = prop_scores$is_m27)
curve_data <- fortify(precrec_obj)
# Plot both Precision-Recall and ROC
plot(precrec_obj)

# Estimate the AUROC and AUPRC in the complete network or the top 200 node subnetwork == M27
prop_scores <- as.data.frame(RWR_M0_Results$RWRM_Results)
prop_scores$Score <- prop_scores$Score / sum(prop_scores$Score)
prop_scores_top200 <- prop_scores[prop_scores$NodeNames %in% top_50_genes, ]

prop_scores_top200$is_m27 <- ifelse(prop_scores_top200$NodeNames %in% M27_genes$V1, 1, 0)
prop_scores_top200$cum_m27 <- cumsum(prop_scores_top200$is_m27)
prop_scores_top200 <- prop_scores_top200[order(prop_scores_top200$Score, decreasing = TRUE), ]
prop_scores_top200$rank <- 1:nrow(prop_scores_top200)

ggplot(prop_scores_top200, aes(x = rank, y = cum_m27)) +
  geom_line(color = "firebrick", size = 1) +
  geom_abline(intercept = 0, slope = sum(prop_scores_top200$is_m27)/nrow(prop_scores_top200), linetype = "dashed") +
  labs(title = "M27 relapse enrichment curve",
       subtitle = "Does the drug perturbation (M0) point toward the relapse (M27) TFs?",
       x = "Ranked genes (by RWR score)",
       y = "Cumulative M27 gemes found") +
  theme_minimal()

precrec_obj <- evalmod(scores = prop_scores_top200$Score, labels = prop_scores_top200$is_m27)

# Plot both Precision-Recall and ROC
plot(precrec_obj)

# Estimate the AUROC and AUPRC in the complete network or the top 200 node subnetwork == M3

prop_scores <- as.data.frame(RWR_M0_Results$RWRM_Results)
prop_scores$Score <- prop_scores$Score / sum(prop_scores$Score)
prop_scores_top200 <- prop_scores[prop_scores$NodeNames %in% top_50_genes, ]

prop_scores_top200$is_m3 <- ifelse(prop_scores_top200$NodeNames %in% M3_genes$V1, 1, 0)
prop_scores_top200$cum_m3 <- cumsum(prop_scores_top200$is_m3)
prop_scores_top200 <- prop_scores_top200[order(prop_scores_top200$Score, decreasing = TRUE), ]
prop_scores_top200$rank <- 1:nrow(prop_scores_top200)

ggplot(prop_scores_top200, aes(x = rank, y = cum_m3)) +
  geom_line(color = "firebrick", size = 1) +
  geom_abline(intercept = 0, slope = sum(prop_scores_top200$is_m27)/nrow(prop_scores_top200), linetype = "dashed") +
  labs(title = "M3 response enrichment curve",
       subtitle = "Does the drug perturbation (M0) point toward the response (M3) genes?",
       x = "Ranked genes (by RWR score)",
       y = "Cumulative M3 gemes found") +
  theme_minimal()

precrec_obj <- evalmod(scores = prop_scores_top200$Score, labels = prop_scores_top200$is_m3)

# Plot both Precision-Recall and ROC
plot(precrec_obj)

### Same deeper analysis 
scores <- prop_scores$Score
labels <- prop_scores$is_m27

precrec_obj <- evalmod(scores = scores, labels = labels)
auc_vals <- auc(precrec_obj)

P <- sum(labels == 1)
N <- sum(labels == 0)
baseline_precision <- P / (P + N)

cat("Positives:", P, "\n")
cat("Negatives:", N, "\n")
cat("Baseline precision:", baseline_precision, "\n")

curve_data <- fortify(precrec_obj)
pr_data <- subset(curve_data, curvetype == "PRC")

max_precision <- max(pr_data$y, na.rm = TRUE)
best_prec_row <- pr_data[which.max(pr_data$y), ]

cat("Max precision:", max_precision, "\n")
cat("Recall at max precision:", best_prec_row$x, "\n")

pr_data$F1 <- 2 * (pr_data$x * pr_data$y) / 
  (pr_data$x + pr_data$y)

max_f1 <- max(pr_data$F1, na.rm = TRUE)
best_f1_row <- pr_data[which.max(pr_data$F1), ]

cat("Max F1:", max_f1, "\n")
cat("Recall at max F1:", best_f1_row$x, "\n")
cat("Precision at max F1:", best_f1_row$y, "\n")

df_ranked <- data.frame(score = scores,
                        label = labels) %>%
  arrange(desc(score)) %>%
  mutate(rank = row_number(),
         cum_tp = cumsum(label),
         precision_at_k = cum_tp / rank)
top_k_values <- c(10, 25, 50, 100, 200)

for (k in top_k_values) {
  prec_k <- df_ranked$precision_at_k[k]
  cat("Precision@", k, ":", prec_k, "\n")
}

for (k in top_k_values) {
  prec_k <- df_ranked$precision_at_k[k]
  enrichment <- prec_k / baseline_precision
  cat("Fold enrichment@", k, ":", enrichment, "\n")
}

df_ranked$recall_at_k <- df_ranked$cum_tp / P

pr_plot <- ggplot(df_ranked, aes(x = rank, y = recall_at_k)) +
  geom_line() +
  labs(title = "Cumulative Recall vs Rank",
       x = "Rank",
       y = "Recall") +
  theme_classic()

pr_plot <- ggplot(pr_data, aes(x = x, y = y)) +
  geom_line(size = 1, color = "darkgoldenrod4") +
  geom_hline(yintercept = baseline_precision,
             linetype = "dashed",
             color = "red") +
  geom_point(aes(x = best_prec_row$x,
                 y = best_prec_row$y),
             color = "blue",
             size = 3) +
  annotate("text",
           x = best_prec_row$x,
           y = best_prec_row$y,
           label = paste0("Max Prec = ",
                          round(max_precision,3)),
           hjust = -0.1) +
  labs(title = "Precision-Recall Curve",
       x = "Recall",
       y = "Precision") +
  theme_classic()
ggsave(pr_plot, filename = "Plots/CLL_only/RWR_from_M0_to_M27.svg", device = "svg", width = 4, height = 4, dpi = 200)

