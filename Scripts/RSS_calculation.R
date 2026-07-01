
library(dplyr)
library(tidyverse)
library(ggplot2)
library(SCENIC)
library(ggrepel)

### ------------------------------------------------------------
### Load the AUCell object and annotate the cells
### ------------------------------------------------------------
AUCell_df <- read.csv("Results/CLL_only/AUcell_CLL_only.csv", header = T, row.names = 1, sep = "\t")

colnames(AUCell_df) <- gsub("\\.\\.\\.", "", colnames(AUCell_df))

cell_names <- rownames(AUCell_df)

M_values <- sub(".*_M([0-9]+)_.*", "\\1", cell_names)
D_values <- paste0("M", M_values)
names(D_values) <- rownames(AUCell_df)
annotation_col <- data.frame(TimePoint = D_values)
rownames(annotation_col) <- cell_names

### ------------------------------------------------------------
### Calculate the RSS scores
### ------------------------------------------------------------
rss <- calcRSS(t(AUCell_df), t(D_values))

rss_df <- rss_2groups %>%
  as.data.frame()
colnames(rss_df) <- rownames(AUCell_df)
rss_df <- rss_df %>% tibble::rownames_to_column("Regulon") %>%
  tidyr::pivot_longer(
    cols = -Regulon,
    names_to = "Cell",
    values_to = "RSS"
  ) %>%
  mutate(Group = D_values[Cell]) %>%
  group_by(Regulon, Group) %>%
  summarise(
    RSS = mean(RSS, na.rm = TRUE),
    .groups = "drop"
  )

rss_df$Group <- factor(rss_df$Group, levels = c("M0", "M3", "M27"))

rss_df <- data.frame(
  regulon = rownames(rss),
  RSS = rss[,1]
)

### ------------------------------------------------------------
### Rank regulons by RSS
### ------------------------------------------------------------

rss_df <- rss_df %>% arrange(desc(RSS))

# Select top regulons to label
topN <- 20
rss_2groups_df$label <- ifelse(1:nrow(rss_2groups_df) <= topN, rss_2groups_df$regulon, "")

tfs_to_label <- up_TFs_M27
rss_df$label <- ifelse(rss_df$regulon %in% label_TFs,
                       rss_df$regulon,
                       "")

### ------------------------------------------------------------
### Plot RSS ranking curve
### ------------------------------------------------------------
rss_df <- data.frame(
  TF  = names(rss_vec),
  RSS = as.numeric(rss_vec)
) %>%
  arrange(desc(RSS)) %>%
  mutate(rank = row_number())
ggplot(rss_df, aes(x = rank, y = RSS)) +
  geom_point(color = "#2C7BB6", size = 2) +
  geom_point(
    data = subset(rss_df, TF %in% tfs_to_label),
    color = "red",
    size = 2.5
  ) +
  geom_text_repel(
    data = subset(rss_df, TF %in% tfs_to_label),
    aes(label = TF),
    color = "red",
    size = 4,
    box.padding = 0.4,
    max.overlaps = Inf
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "Regulon",
    y = "RSS"
  )


RSS_plot <- ggplot(rss_df, aes(x = seq_along(RSS), y = RSS)) +
  geom_point(color = "steelblue", size = 2) +
  geom_point(data = rss_df[1:topN,], aes(x = seq_len(topN), y = RSS),
             color = "red", size = 3) +
  geom_text_repel(
    #data = rss_df[1:topN,],
    data = rss_df %>% filter(regulon %in% label_TFs),
    #aes(label = label, x = seq_len(topN), y = RSS),
    aes(x = seq_along(RSS)[regulon %in% label_TFs], y = RSS[regulon %in% label_TFs], label = regulon),
    color = "red",
    nudge_x = 0.5,
    nudge_y = 0.01,
    max.overlaps = Inf,
    box.padding = 0.4,
    force = 2,
    segment.color = "grey50",
    direction = "y",
    seed = 42
  ) +
  labs(x = "Regulon", y = "RSS", title = "RSS score, M27") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank()
  )

ggsave(filename = "../Plots_Results/RSS_plot_M27.svg", plot = RSS_plot, device = "svg", width = 5, height = 7)

### ------------------------------------------------------------
### Violin plots of top 5 TF regulons
### ------------------------------------------------------------
n_TFs_of_interest = 9

top_TFs <- rss_M27_df %>% 
  arrange(desc(RSS)) %>% 
  slice_head(n = n_TFs_of_interest) %>%
  pull(Regulon)

mat_top5 <- AUCell_df[, top_TFs, drop = FALSE]

df_plot <- mat_top5 %>%
  as.data.frame() %>%
  mutate(Cell = rownames(.),
         Group = M_values[Cell]) %>%   # add M0/M3/M27
  pivot_longer(cols = top_TFs,
               names_to = "Regulon",
               values_to = "Activity")
df_plot$Group <- factor(groups, levels = c("M0", "M3", "M27"))

plot_regulon_distribution <- ggplot(df_plot, aes(x = Group, y = Activity, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.5) +
  facet_wrap(~ Regulon, scales = "free_y", nrow = 4) +
  theme_bw(base_size = 14) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Cell category", y = "Regulon activity",
       title = paste0("Top ", n_TFs_of_interest, " regulons in M27 category")) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 12, face = "bold")
  )

ggsave("Plots/CLL_only/top_9_regulon_AUCell.svg", plot = plot_regulon_distribution, width = 8, height = 6, dpi = 300)

### ------------------------------------------------------------
### Perform Wilcoxon rank-sum test to identify TF with a high regulon activity in M27
### ------------------------------------------------------------

library(tidyr)

df_long <- AUCell_df %>%
  as.data.frame() %>%
  mutate(Cell = rownames(.),
         Group = D_values[Cell]) %>%
  pivot_longer(
    cols = everything()[!colnames(.) %in% c("Cell","Group")],
    names_to = "Regulon",
    values_to = "Activity"
  )

df_long$Group <- factor(df_long$Group, levels = c("M0","M3","M27"))

results <- df_long %>%
  group_by(Regulon) %>%
  summarise(
    p_value = wilcox.test(
      Activity[Group == "M27"],
      Activity[Group %in% c("M0","M3")]
    )$p.value,
    direction = median(Activity[Group == "M27"]) -
      median(Activity[Group %in% c("M0","M3")])
  ) %>%
  mutate(
    FDR = p.adjust(p_value, method = "fdr")
  ) %>%
  arrange(FDR)

sig_up_M27   <- results %>% filter(FDR < 0.05 & direction > 0 & p_value < 0.05)
sig_down_M27 <- results %>% filter(FDR < 0.05 & direction < 0 & p_value < 0.05)

new_labels <- ifelse(D_values == "M27","M27","Other")
names(new_labels) <- names(D_values)

rss_2groups <- calcRSS(t(AUCell_df), t(new_labels))
