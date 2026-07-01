###############################################
###   CENTRALITY ANALYSIS FOR GRN NETWORK   ###
###############################################

# Install packages if needed
packages <- c("igraph", "tidyverse")
installed <- packages %in% rownames(installed.packages())
if (any(!installed)) install.packages(packages[!installed])

library(igraph)
library(tidyverse)
library(patchwork)
library(viridis)

###############################################
### 1. Load the GRN Table (Source–Target)   ###
###############################################

# Path to your GRN file
grn_file <- "Results/CLL_only/SCENIC_GRN_edges.csv"

# Read network (assumed columns: Source, Target, Mode)
grn <- read.csv(grn_file, stringsAsFactors = FALSE)

# Inspect
head(grn)
unique(grn$Mode)  # e.g., activation, inhibition, unknown

###############################################
### 2. Build a Directed Graph in igraph     ###
###############################################

# Build graph from the Source and Target columns
g <- graph_from_data_frame(
  d = grn[, c("Source", "Target")], 
  directed = TRUE
)

# Optional: add Mode as edge attribute
E(g)$mode <- grn$Mode

###############################################
### 3. Compute Centrality Measures          ###
###############################################

# Degree centrality
degree_total  <- degree(g, mode = "all")
degree_in     <- degree(g, mode = "in")
degree_out    <- degree(g, mode = "out")

# Betweenness centrality
betw <- betweenness(g, directed = TRUE)

# Closeness centrality
close <- closeness(g, mode = "out")

# Eigenvector centrality
eig <- eigen_centrality(g, directed = TRUE, scale = TRUE)$vector

# PageRank
pr <- page.rank(g, directed = TRUE)$vector

###############################################
### 4. Combine All Centrality Results       ###
###############################################

centrality_df <- data.frame(
  Gene = names(degree_total),
  Degree = degree_total,
  In_Degree = degree_in,
  Out_Degree = degree_out,
  Betweenness = betw,
  Closeness = close,
  Eigenvector = eig,
  PageRank = pr
)

# Rank nodes by each centrality measure
centrality_df_ranked <- centrality_df %>%
  arrange(desc(Degree))

###############################################
### 5. Display Top Central Genes            ###
###############################################

top_n <- 20  # how many top genes to show

cat("\n===== TOP GENES BY DEGREE =====\n")
print(centrality_df %>% arrange(desc(Degree)) %>% head(top_n))

cat("\n===== TOP GENES BY BETWEENNESS =====\n")
print(centrality_df %>% arrange(desc(Betweenness)) %>% head(top_n))

cat("\n===== TOP GENES BY PAGERANK =====\n")
print(centrality_df %>% arrange(desc(PageRank)) %>% head(top_n))

###############################################
### 6. Plotting (Optional)                  ###
###############################################

# Histogram of degree distribution
hist(degree_total, 
     breaks = 30, 
     main = "Degree Distribution",
     xlab = "Degree")

theme_pub <- theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.6),
    plot.margin = margin(10, 10, 10, 10)
  )

# Top Genes based on a centrality measure
p_deg <- centrality_df %>%
  arrange(desc(Degree)) %>% slice(1:top_n) %>%
  ggplot(aes(x = reorder(Gene, Degree), y = Degree, fill = Degree)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  coord_flip() +
  scale_fill_viridis(option = "D", direction = 1) +
  labs(title = "Top 20 Genes by Degree") +
  theme_pub +
  theme(legend.position = "none")

### Top 20 by Betweenness
p_betw <- centrality_df %>%
  arrange(desc(Betweenness)) %>% slice(1:top_n) %>%
  ggplot(aes(x = reorder(Gene, Betweenness), y = Betweenness, fill = Betweenness)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  coord_flip() +
  scale_fill_viridis(option = "C", direction = 1) +
  labs(title = "Top 20 Genes by Betweenness") +
  theme_pub +
  theme(legend.position = "none")

### Top 20 by PageRank
p_pr <- centrality_df %>%
  arrange(desc(PageRank)) %>% slice(1:top_n) %>%
  ggplot(aes(x = reorder(Gene, PageRank), y = PageRank, fill = PageRank)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  coord_flip() +
  scale_fill_viridis(option = "A", direction = 1) +
  labs(title = "Top 20 Genes by PageRank") +
  theme_pub +
  theme(legend.position = "none")

combined_plot <- p_deg | p_betw
combined_plot

ggsave("Plots/CLL_only/top_degree_TFs.svg", plot = combined_plot, width = 9, height = 6, dpi = 300)

###############################################
### 7. Save Outputs                         ###
###############################################

write.csv(centrality_df, "GRN_centrality_measures.csv", row.names = FALSE)
write.csv(centrality_df_ranked, "GRN_centrality_ranked.csv", row.names = FALSE)

