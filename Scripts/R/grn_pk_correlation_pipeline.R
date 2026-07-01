## ============================================================
## GRN-PK ibrutinib resistance study — full pipeline
## Approach: PK-gene expression Spearman correlation
##   (Option 5: stratify genes by PK sensitivity trajectory)
##
## PK data: ID6_AUC_comparaison_methodes_final.xlsx
##   - day           : absolute time in days (1 row per day)
##   - Cmax_IBRU     : daily predicted Cmax (ng/mL)
##   - AUC_IBRU_luld_pknca : daily AUC (linear-up log-down, PKNCA)
##   - AUC_cum_IBRU  : cumulative AUC from day 1 (ng/mL·day)
##   - imputed       : FALSE = observed month, TRUE = imputed
##
## Key timepoints:
##   M3  = month 3 (day ~90)  — response
##   M27 = month 27 (day ~810) — relapse
##
## Outputs: all figures in SVG + PDF, all tables as CSV
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(readxl)
  library(tidyverse)
  library(ggrepel)
  library(patchwork)
  library(pheatmap)
  library(viridis)
  library(caret)
  library(h2o)
  library(svglite)   # install.packages("svglite") if missing
})

## ============================================================
## 0. CONFIGURATION
## ============================================================

SEURAT_PATH     <- "Data/CLL_only/cll_obj_tfAct_AUCell.rds"
PK_FILE         <- "Results/CLL_only/PK_results/ID6_AUC_comparaison_methodes_final.xlsx"
PK_SHEET        <- 1
DRIVER_M3_PATH  <- "Results/CLL_only/CellRank_results/driver_genes_cellfate3.csv"
DRIVER_M27_PATH <- "Results/CLL_only/CellRank_results/driver_genes_cellfate4.csv"
OUTPUT_DIR      <- "Results/CLL_only/ML_results/results_pk_correlation"

# Clinical timepoints
M3_DAY          <- 90
M27_DAY         <- 822

# AutoML settings (used for the baseline R² validation block)
N_REPLICATES    <- 10
TEST_FRACTION   <- 0.20
AUTOML_TIME_SEC <- 60
H2O_MAX_MEM     <- "4G"
CV_FOLDS        <- 5
SEED            <- 42
NOISE_CV        <- 0.10

# Figure font size (increased for publication)
BASE_FONT       <- 14

dir.create(OUTPUT_DIR, showWarnings = FALSE)

## ============================================================
## HELPER: save figure in both SVG and PDF
## ============================================================

save_fig <- function(plot, filename, width = 10, height = 7) {
  base <- file.path(OUTPUT_DIR, filename)
  ggsave(paste0(base, ".svg"), plot,
         width = width, height = height, device = svglite)
  ggsave(paste0(base, ".pdf"), plot,
         width = width, height = height, device = cairo_pdf)
  message(sprintf("  Saved: %s.svg / .pdf", base))
}

theme_pub <- theme_minimal(base_size = BASE_FONT) +
  theme(
    plot.title       = element_text(face = "bold", size = BASE_FONT + 2),
    plot.subtitle    = element_text(size = BASE_FONT - 1, colour = "grey40"),
    strip.text       = element_text(face = "bold", size = BASE_FONT),
    legend.position  = "bottom",
    legend.text      = element_text(size = BASE_FONT - 1),
    legend.title     = element_text(size = BASE_FONT - 1, face = "bold"),
    axis.title       = element_text(size = BASE_FONT),
    axis.text        = element_text(size = BASE_FONT - 1),
    panel.grid.minor = element_blank()
  )

class_cols <- c(
  "Acquired resistance" = "#d95f02",
  "Intrinsic"          = "#1f78b4",
  "Response-specific"  = "#33a02c",
  "Neutral"            = "grey78",
  "PK-sensitive"       = "#9932CC",
  "PK-insensitive"     = "grey78"
)

## ============================================================
## 1. LOAD PK DATA
## ============================================================

message("\n[1/8] Loading PK data...")

pk_raw <- read_excel(PK_FILE, sheet = PK_SHEET)

pk_df <- pk_raw %>%
  filter(!is.na(day), day >= 0, day <= M27_DAY)

# Separate observed vs imputed for plotting
pk_obs  <- pk_df %>% filter(imputed == FALSE)
pk_all  <- pk_df   # use all rows (obs + imputed) for cumulative AUROC

message(sprintf("  Total rows: %d  |  Observed: %d  |  Imputed: %d",
                nrow(pk_df), nrow(pk_obs), nrow(pk_df) - nrow(pk_obs)))

## ============================================================
## 2. EXTRACT PK SUMMARY VALUES AT M3 AND M27
## ============================================================

message("\n[2/8] Extracting PK values at M3 and M27...")

get_pk_at <- function(df, target_day, col) {
  df %>%
    mutate(dist = abs(day - target_day)) %>%
    slice_min(dist, n = 1, with_ties = FALSE) %>%
    pull({{ col }})
}

# Cmax at M3 and M27
Cmax_M3  <- get_pk_at(pk_all, M3_DAY,  Cmax_IBRU)
Cmax_M27 <- get_pk_at(pk_all, M27_DAY, Cmax_IBRU)

# Daily AUC (luld pknca) at M3 and M27
AUC_M3   <- get_pk_at(pk_all, M3_DAY,  AUC_IBRU_luld_pknca)
AUC_M27  <- get_pk_at(pk_all, M27_DAY, AUC_IBRU_luld_pknca)

# Cumulative AUC at M3 and M27
AUCcum_M3  <- get_pk_at(pk_all, M3_DAY,  AUC_cum_IBRU)
AUCcum_M27 <- get_pk_at(pk_all, M27_DAY, AUC_cum_IBRU)

# PK variability (CV%) per window
pk_cv_M3 <- pk_all %>%
  filter(day <= M3_DAY) %>%
  summarise(cv = sd(Cmax_IBRU, na.rm = TRUE) /
                 mean(Cmax_IBRU, na.rm = TRUE) * 100) %>%
  pull(cv)

pk_cv_M27 <- pk_all %>%
  filter(day <= M27_DAY) %>%
  summarise(cv = sd(Cmax_IBRU, na.rm = TRUE) /
                 mean(Cmax_IBRU, na.rm = TRUE) * 100) %>%
  pull(cv)

pk_summary <- tibble(
  metric   = c("Cmax (ng/mL)", "Daily AUC (ng/mL·day)",
               "Cumulative AUC (ng/mL·day)", "CV% Cmax"),
  M3       = c(Cmax_M3,  AUC_M3,  AUCcum_M3,  pk_cv_M3),
  M27      = c(Cmax_M27, AUC_M27, AUCcum_M27, pk_cv_M27)
)

cat("\n===== PK SUMMARY VALUES =====\n")
print(pk_summary)
write_csv(pk_summary, file.path(OUTPUT_DIR, "pk_summary_values.csv"))

## ============================================================
## 3. PK PROFILE FIGURE
## ============================================================

message("\n[3/8] Generating PK profile figure...")

# Daily mean Cmax for simulated curve (all days)
pk_daily_pred <- pk_all %>%
  mutate(day_round = round(day)) %>%
  group_by(day_round) %>%
  summarise(
    mean_Cmax    = mean(Cmax_IBRU,          na.rm = TRUE),
    mean_AUCcum  = mean(AUC_cum_IBRU,       na.rm = TRUE),
    .groups = "drop"
  )

# Observed: daily stats (observed months only)
pk_daily_obs <- pk_obs %>%
  mutate(day_round = floor(day)) %>%
  group_by(day_round) %>%
  summarise(
    median_Cmax = median(Cmax_IBRU, na.rm = TRUE),
    q25_Cmax    = quantile(Cmax_IBRU, 0.25, na.rm = TRUE),
    q75_Cmax    = quantile(Cmax_IBRU, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

y_max <- max(pk_daily_obs$q75_Cmax, na.rm = TRUE) * 1.25

p_pk <- ggplot() +
  # Cumulative AUROC shading
  annotate("rect", xmin = 0, xmax = M27_DAY, ymin = 0, ymax = Inf,
           fill = "#d95f02", alpha = 0.05) +
  annotate("rect", xmin = 0, xmax = M3_DAY, ymin = 0, ymax = Inf,
           fill = "#2171b5", alpha = 0.10) +
  # Observed IQR ribbon and median (observed months only)
  geom_ribbon(data = pk_daily_obs,
              aes(x = day_round, ymin = q25_Cmax, ymax = q75_Cmax),
              fill = "#2171b5", alpha = 0.22) +
  geom_point(data = pk_daily_obs,
             aes(x = day_round, y = median_Cmax),
             colour = "#2171b5", size = 1.2, alpha = 0.80) +
  # Last observed line
  geom_vline(xintercept = max(pk_obs$day, na.rm = TRUE),
             linetype = "dotted", colour = "#2171b5", linewidth = 0.8) +
  annotate("text",
           x = max(pk_obs$day, na.rm = TRUE) + 10,
           y = y_max * 0.70,
           label = "Last observed\nconcentration (M6)",
           hjust = 0, size = BASE_FONT / 3.5, colour = "#2171b5") +
  # Simulated daily mean curve (full window)
  geom_line(data = pk_daily_pred,
            aes(x = day_round, y = mean_Cmax),
            colour = "#d95f02", linewidth = 0.9, alpha = 0.90) +
  # Timepoint lines
  geom_vline(xintercept = c(0, M3_DAY, M27_DAY),
             linetype = "dashed", colour = "grey35", linewidth = 0.6) +
  annotate("text",
           x = c(2, M3_DAY + 3, M27_DAY - 10),
           y = y_max,
           vjust = 1.6, hjust = 0,
           label = c("M0", "M3", "M27"),
           size = BASE_FONT / 3, colour = "grey25", fontface = "italic") +
  # AUROC annotations
  annotate("text", x = M3_DAY / 2, y = y_max * 0.92,
           label = sprintf("Cumul. AUROC\n0\u2192M3\n%.0f ng/mL\u00b7day", AUCcum_M3),
           size = BASE_FONT / 3.5, colour = "#2171b5",
           fontface = "bold", hjust = 0.5) +
  annotate("text", x = (M3_DAY + M27_DAY) / 2, y = y_max * 0.92,
           label = sprintf("Cumul. AUROC\n0\u2192M27\n%.0f ng/mL\u00b7day", AUCcum_M27),
           size = BASE_FONT / 3.5, colour = "#d95f02",
           fontface = "bold", hjust = 0.5) +
  coord_cartesian(ylim = c(0, y_max)) +
  scale_x_continuous(
    name   = "Time (days)",
    breaks = c(0, 90, 180, 270, 360, 450, 540, 630, 720, 810),
    labels = c("0\n(M0)", "90\n(M3)", "180\n(M6)",
               "270", "360", "450", "540", "630", "720", "810\n(M27)")
  ) +
  scale_y_continuous(name = "Ibrutinib Cmax (ng/mL)") +
  labs(
    title    = "Ibrutinib PK profile — observed vs simulated",
    subtitle = paste0(
      "Blue: observed daily median \u00b1 IQR (M1\u2013M6).  ",
      "Red: simulated daily mean Cmax (M0\u2013M27).\n",
      "Shaded: cumulative AUROC windows (blue = 0\u2192M3; orange = 0\u2192M27)."
    )
  ) +
  theme_pub

save_fig(p_pk, "fig1_pk_profile", width = 12, height = 6)

## ============================================================
## 4. LOAD SEURAT + AUCCELL
## ============================================================

message("\n[4/8] Loading Seurat object and AUCell data...")

seu <- readRDS(SEURAT_PATH)
DefaultAssay(seu) <- "RNA"

get_pseudobulk <- function(seurat_obj, tp) {
  cells <- colnames(seurat_obj)[seurat_obj$time_point == tp]
  rowMeans(GetAssayData(seurat_obj, layer = "data")[, cells, drop = FALSE])
}
get_pseudobulk_sd <- function(seurat_obj, tp) {
  cells <- colnames(seurat_obj)[seurat_obj$time_point == tp]
  apply(GetAssayData(seurat_obj, layer = "data")[, cells, drop = FALSE],
        1, sd, na.rm = TRUE)
}

expr_M0  <- get_pseudobulk(seu, "M0")
expr_M3  <- get_pseudobulk(seu, "M3")
expr_M27 <- get_pseudobulk(seu, "M27")
sd_M0    <- get_pseudobulk_sd(seu, "M0")

# AUCell from Seurat assay (TF × cell matrix)
DefaultAssay(seu) <- "AUCell"
auc_data <- GetAssayData(seu, layer = "data")

get_auc_pseudobulk <- function(mat, tp) {
  cells <- intersect(colnames(seu)[seu$time_point == tp], colnames(mat))
  rowMeans(mat[, cells, drop = FALSE])
}
get_auc_sd <- function(mat, tp) {
  cells <- intersect(colnames(seu)[seu$time_point == tp], colnames(mat))
  apply(mat[, cells, drop = FALSE], 1, sd, na.rm = TRUE)
}

auc_M0    <- get_auc_pseudobulk(auc_data, "M0")
auc_M3    <- get_auc_pseudobulk(auc_data, "M3")
auc_M27   <- get_auc_pseudobulk(auc_data, "M27")
auc_sd_M0 <- get_auc_sd(auc_data, "M0")
tfs_use   <- rownames(auc_data)

DefaultAssay(seu) <- "RNA"
message(sprintf("  TFs: %d  |  Genes: %d", length(tfs_use), length(expr_M0)))

## ============================================================
## 5. DRIVER GENES
## ============================================================

drivers_M3  <- read_csv(DRIVER_M3_PATH) %>% as.data.frame()
drivers_M3 <- drivers_M3 %>% filter(drivers_M3$`3_corr`>0.1)
drivers_M3 <- drivers_M3[["...1"]]

drivers_M27 <- read_csv(DRIVER_M27_PATH) %>% as.data.frame()
drivers_M27 <- drivers_M27 %>% filter(drivers_M27$`4_corr`>0.1)
drivers_M27 <- drivers_M27[["...1"]]

genes_use   <- intersect(union(drivers_M3, drivers_M27), rownames(seu))
message(sprintf("  Driver genes: %d", length(genes_use)))

## ============================================================
## 6. OPTION 5 — PK-GENE EXPRESSION CORRELATION
##
## For each gene and TF, compute the Spearman correlation between
## its expression trajectory (M0, M3, M27) and the cumulative
## AUROC trajectory (0, AUCcum_M3, AUCcum_M27).
## Genes with high |rho| track drug exposure — PK-sensitive.
## Genes with low |rho| are PK-independent — intrinsic programme.
## ============================================================

message("\n[5/8] Computing PK-expression trajectory correlations...")

pk_traj <- c(0, AUCcum_M3, AUCcum_M27)   # cumulative AUROC at M0, M3, M27
message(sprintf("  PK trajectory (cumul. AUROC): M0=0  M3=%.0f  M27=%.0f ng/mL·day",
                AUCcum_M3, AUCcum_M27))

# Also compute with Cmax trajectory
cmax_traj <- c(0, Cmax_M3, Cmax_M27)
message(sprintf("  Cmax trajectory:              M0=0  M3=%.1f  M27=%.1f ng/mL",
                Cmax_M3, Cmax_M27))

# -- Gene correlation
gene_pk_cor <- tibble(
  gene     = genes_use,
  expr_M0  = expr_M0[genes_use],
  expr_M3  = expr_M3[genes_use],
  expr_M27 = expr_M27[genes_use],
  sd_M0    = sd_M0[genes_use]
) %>%
  filter(expr_M0 > 0) %>%
  rowwise() %>%
  mutate(
    # Spearman correlation: expression trajectory vs cumulative AUROC
    rho_AUCcum = cor(
      c(expr_M0, expr_M3, expr_M27), pk_traj,
      method = "spearman"
    ),
    # Spearman correlation: expression trajectory vs Cmax
    rho_Cmax = cor(
      c(expr_M0, expr_M3, expr_M27), cmax_traj,
      method = "spearman"
    ),
    # Expression change magnitude (for ranking)
    delta_M3  = expr_M3  - expr_M0,
    delta_M27 = expr_M27 - expr_M0,
    total_change = abs(delta_M3) + abs(delta_M27)
  ) %>%
  ungroup() %>%
  # PK sensitivity classification
  mutate(
    pk_class = case_when(
      abs(rho_AUCcum) >= 0.9 & rho_AUCcum  > 0 ~ "PK-driven (activated)",
      abs(rho_AUCcum) >= 0.9 & rho_AUCcum  < 0 ~ "PK-driven (suppressed)",
      abs(rho_AUCcum) <  0.5                    ~ "PK-independent",
      TRUE                                       ~ "Partial PK correlation"
    )
  ) %>%
  arrange(desc(abs(rho_AUCcum)))

# -- TF correlation
tf_pk_cor <- tibble(
  tf       = tfs_use,
  auc_M0   = auc_M0[tfs_use],
  auc_M3   = auc_M3[tfs_use],
  auc_M27  = auc_M27[tfs_use],
  sd_M0    = auc_sd_M0[tfs_use]
) %>%
  filter(auc_M0 > 0) %>%
  rowwise() %>%
  mutate(
    rho_AUCcum = cor(
      c(auc_M0, auc_M3, auc_M27), pk_traj,
      method = "spearman"
    ),
    rho_Cmax = cor(
      c(auc_M0, auc_M3, auc_M27), cmax_traj,
      method = "spearman"
    ),
    delta_M3  = auc_M3  - auc_M0,
    delta_M27 = auc_M27 - auc_M0,
    total_change = abs(delta_M3) + abs(delta_M27)
  ) %>%
  ungroup() %>%
  mutate(
    pk_class = case_when(
      abs(rho_AUCcum) >= 0.9 & rho_AUCcum  > 0 ~ "PK-driven (activated)",
      abs(rho_AUCcum) >= 0.9 & rho_AUCcum  < 0 ~ "PK-driven (suppressed)",
      abs(rho_AUCcum) <  0.5                    ~ "PK-independent",
      TRUE                                       ~ "Partial PK correlation"
    )
  ) %>%
  arrange(desc(abs(rho_AUCcum)))

write_csv(gene_pk_cor, file.path(OUTPUT_DIR, "gene_pk_correlation.csv"))
write_csv(tf_pk_cor,   file.path(OUTPUT_DIR, "tf_pk_correlation.csv"))

cat("\n--- Gene PK class distribution ---\n")
print(table(gene_pk_cor$pk_class))
cat("\n--- TF PK class distribution ---\n")
print(table(tf_pk_cor$pk_class))

cat("\n--- Top 15 PK-sensitive genes (|rho| by cumulative AUROC) ---\n")
print(gene_pk_cor %>% select(gene, rho_AUCcum, rho_Cmax, pk_class) %>% head(15))

cat("\n--- Top 15 PK-sensitive TFs ---\n")
print(tf_pk_cor %>% select(tf, rho_AUCcum, rho_Cmax, pk_class) %>% head(15))

## ============================================================
## 7. AUTOML — BASELINE R² VALIDATION
##    Run AutoML with the real PK values from the Excel file.
##    This produces the ΔR² comparison figure for the manuscript.
## ============================================================

message("\n[6/8] h2o AutoML — baseline R² validation...")

h2o.init(max_mem_size = H2O_MAX_MEM, nthreads = -1)
h2o.no_progress()

# Feature matrix builder (uses real PK values from Excel)
make_gene_X <- function(df, task) {
  pk_val <- if (task == "M3") AUCcum_M3 else AUCcum_M27
  cmax   <- if (task == "M3") Cmax_M3   else Cmax_M27
  n      <- nrow(df)
  set.seed(SEED + n)   # reproducible noise per task
  df %>% transmute(
    expr_M0   = expr_M0,
    sd_M0     = sd_M0,
    AUCcum    = pk_val * (1 + rnorm(n, 0, NOISE_CV)),
    Cmax      = cmax   * (1 + rnorm(n, 0, NOISE_CV)),
    intx_cum  = expr_M0 * AUCcum,
    intx_cmax = expr_M0 * Cmax
  )
}

make_tf_X <- function(df, task) {
  pk_val <- if (task == "M3") AUCcum_M3 else AUCcum_M27
  cmax   <- if (task == "M3") Cmax_M3   else Cmax_M27
  n      <- nrow(df)
  set.seed(SEED + n + 1)
  df %>% transmute(
    auc_M0    = auc_M0,
    auc_sd_M0 = sd_M0,
    AUCcum    = pk_val * (1 + rnorm(n, 0, NOISE_CV)),
    Cmax      = cmax   * (1 + rnorm(n, 0, NOISE_CV)),
    intx_cum  = auc_M0 * AUCcum,
    intx_cmax = auc_M0 * Cmax
  )
}

remove_baseline <- function(train_df, test_df, target, predictor) {
  fit <- lm(as.formula(paste(target, "~", predictor)), data = train_df)
  list(
    train_resid = residuals(fit),
    test_resid  = test_df[[target]] - predict(fit, newdata = test_df),
    baseline_r2 = summary(fit)$r.squared
  )
}

add_noise <- function(df, noise_cv, seed) {
  set.seed(seed)
  df$AUROC_M3  <- AUCcum_M3  * (1 + rnorm(nrow(df), 0, noise_cv))
  df$AUROC_M27 <- AUCcum_M27 * (1 + rnorm(nrow(df), 0, noise_cv))
  df
}

# Build base data frames
base_gene_df <- tibble(
  gene     = genes_use,
  expr_M0  = expr_M0[genes_use],
  expr_M3  = expr_M3[genes_use],
  expr_M27 = expr_M27[genes_use],
  sd_M0    = sd_M0[genes_use]
) %>% filter(expr_M0 > 0) %>%
  mutate(expr_decile = ntile(expr_M0, 10))

base_tf_df <- tibble(
  tf        = tfs_use,
  auc_M0    = auc_M0[tfs_use],
  auc_M3    = auc_M3[tfs_use],
  auc_M27   = auc_M27[tfs_use],
  sd_M0     = auc_sd_M0[tfs_use]
) %>% filter(auc_M0 > 0) %>%
  mutate(auc_decile = ntile(auc_M0, 10))

# Train/test split
set.seed(SEED)
idx_gene <- createDataPartition(base_gene_df$expr_decile,
                                 p = 1 - TEST_FRACTION, list = FALSE)[, 1]
train_gene <- base_gene_df[ idx_gene, ]
test_gene  <- base_gene_df[-idx_gene, ]

set.seed(SEED + 1)
idx_tf <- createDataPartition(base_tf_df$auc_decile,
                               p = 1 - TEST_FRACTION, list = FALSE)[, 1]
train_tf <- base_tf_df[ idx_tf, ]
test_tf  <- base_tf_df[-idx_tf, ]

tasks_cfg <- list(
  list(track="gene", task_id="M3",  target="expr_M3",
       pred="expr_M0", train=train_gene, test=test_gene, feat=make_gene_X),
  list(track="gene", task_id="M27", target="expr_M27",
       pred="expr_M0", train=train_gene, test=test_gene, feat=make_gene_X),
  list(track="tf",   task_id="M3",  target="auc_M3",
       pred="auc_M0",  train=train_tf,   test=test_tf,   feat=make_tf_X),
  list(track="tf",   task_id="M27", target="auc_M27",
       pred="auc_M0",  train=train_tf,   test=test_tf,   feat=make_tf_X)
)

all_metrics  <- list()
baseline_r2  <- list()
final_models <- list()

for (rep in seq_len(N_REPLICATES)) {
  message(sprintf("\n  Replicate %d / %d", rep, N_REPLICATES))
  for (cfg in tasks_cfg) {
    label    <- paste0(cfg$track, "_", cfg$task_id)
    seed_rep <- SEED + rep * 100 +
                which(sapply(tasks_cfg, \(x)
                  x$task_id == cfg$task_id && x$track == cfg$track))

    tr <- add_noise(cfg$train, NOISE_CV, seed_rep)
    te <- add_noise(cfg$test,  NOISE_CV, seed_rep + 500)

    bl <- remove_baseline(tr, te, cfg$target, cfg$pred)
    tr$resid <- bl$train_resid
    te$resid <- bl$test_resid
    if (rep == 1) baseline_r2[[label]] <- bl$baseline_r2

    X_tr <- cfg$feat(tr, cfg$task_id)
    X_te <- cfg$feat(te, cfg$task_id)

    tr_h2o <- as.h2o(bind_cols(X_tr, resid = tr$resid))
    te_h2o <- as.h2o(bind_cols(X_te, resid = te$resid))

    aml <- h2o.automl(
      x = colnames(X_tr), y = "resid",
      training_frame     = tr_h2o,
      nfolds             = CV_FOLDS,
      max_runtime_secs   = AUTOML_TIME_SEC,
      seed               = seed_rep,
      sort_metric        = "RMSE",
      stopping_metric    = "RMSE",
      stopping_rounds    = 3,
      keep_cross_validation_predictions = TRUE
    )

    pred_te <- as.vector(h2o.predict(aml@leader, te_h2o))
    r2      <- if (sd(te$resid) > 0) cor(te$resid, pred_te)^2 else NA_real_
    sp      <- cor(te$resid, pred_te, method = "spearman")
    algo    <- str_extract(aml@leader@model_id,
                           "GBM|XGBoost|DRF|XRT|GLM|DeepLearning|StackedEnsemble")

    all_metrics[[length(all_metrics) + 1]] <- tibble(
      track       = cfg$track, task = cfg$task_id,
      replicate   = rep, best_model = algo,
      R2_test     = r2, Spearman = sp,
      baseline_R2 = baseline_r2[[label]],
      delta_R2    = r2 - baseline_r2[[label]]
    )

    if (rep == N_REPLICATES)
      final_models[[label]] <- list(pred = pred_te, y_te = te$resid,
                                    te_base = cfg$test, bl = bl)

    h2o.rm(tr_h2o); h2o.rm(te_h2o)
    cat(sprintf("    [%s] rep%02d  %-16s  R²=%.3f  ΔR²=%+.4f\n",
                label, rep,
                ifelse(is.na(algo), "unknown", algo),
                r2, r2 - baseline_r2[[label]]))
  }
}

metrics_df <- bind_rows(all_metrics)
summary_df <- metrics_df %>%
  group_by(track, task) %>%
  summarise(
    mean_R2       = mean(R2_test,  na.rm = TRUE),
    sd_R2         = sd(R2_test,    na.rm = TRUE),
    mean_delta_R2 = mean(delta_R2, na.rm = TRUE),
    sd_delta_R2   = sd(delta_R2,   na.rm = TRUE),
    mean_Spearman = mean(Spearman, na.rm = TRUE),
    baseline_R2   = first(baseline_R2),
    .groups = "drop"
  )

write_csv(metrics_df,  file.path(OUTPUT_DIR, "automl_metrics.csv"))
write_csv(summary_df,  file.path(OUTPUT_DIR, "automl_summary.csv"))

h2o.shutdown(prompt = FALSE)

## ============================================================
## 8. FIGURES
## ============================================================

message("\n[7/8] Generating figures...")

# -- F2: PK-gene correlation scatter (cumulative AUROC rho)
top_genes_pk <- gene_pk_cor %>%
  filter(pk_class != "PK-independent") %>%
  slice_max(abs(rho_AUCcum), n = 15, with_ties = FALSE)

pk_class_cols <- c(
  "PK-driven (activated)"   = "#d95f02",
  "PK-driven (suppressed)"  = "#1f78b4",
  "Partial PK correlation"  = "#33a02c",
  "PK-independent"          = "grey78"
)

p_gene_cor <- gene_pk_cor %>%
  ggplot(aes(x = rho_AUCcum, y = total_change, colour = pk_class)) +
  geom_point(alpha = 0.60, size = 2.0) +
  geom_vline(xintercept = c(-0.9, 0.9),
             linetype = "dashed", colour = "grey40", linewidth = 0.6) +
  geom_vline(xintercept = c(-0.5, 0.5),
             linetype = "dotted", colour = "grey60", linewidth = 0.5) +
  geom_text_repel(data = top_genes_pk,
                  aes(label = gene), size = BASE_FONT / 3.5,
                  max.overlaps = 20, show.legend = FALSE) +
  scale_colour_manual(values = pk_class_cols, name = "PK class") +
  labs(
    title    = "Gene expression trajectory vs cumulative AUROC",
    #subtitle = "Spearman \u03c1 between expression (M0\u2192M3\u2192M27) and cumulative ibrutinib AUROC",
    x        = "Spearman \u03c1 (expression ~ cumul. AUROC)",
    y        = "Total expression change |ΔM3| + |ΔM27|"
  ) +
  theme_pub + 
  theme(
    axis.title       = element_text(size = 14),
    axis.text        = element_text(size = 12),
    legend.text      = element_text(size = 16),
    plot.title       = element_text(size = 14, face = "bold"),
    legend.title    = element_blank(),
    legend.position  = "bottom"
  )

save_fig(p_gene_cor, "fig2_gene_pk_correlation", width = 8.3, height = 7)

# -- F3: TF PK correlation scatter
top_tfs_pk <- tf_pk_cor %>%
  filter(pk_class != "PK-independent") %>%
  slice_max(abs(rho_AUCcum), n = 15, with_ties = FALSE)

p_tf_cor <- tf_pk_cor %>%
  ggplot(aes(x = rho_AUCcum, y = total_change, colour = pk_class)) +
  geom_point(alpha = 0.60, size = 2.0) +
  geom_vline(xintercept = c(-0.9, 0.9),
             linetype = "dashed", colour = "grey40", linewidth = 0.6) +
  geom_vline(xintercept = c(-0.5, 0.5),
             linetype = "dotted", colour = "grey60", linewidth = 0.5) +
  geom_text_repel(data = top_tfs_pk,
                  aes(label = tf), size = BASE_FONT / 3.5,
                  max.overlaps = 20, show.legend = FALSE) +
  scale_colour_manual(values = pk_class_cols, name = "PK class") +
  labs(
    title    = "TF activity trajectory vs cumulative AUROC",
    subtitle = "Spearman \u03c1 between AUCell score (M0\u2192M3\u2192M27) and cumulative ibrutinib AUROC",
    x        = "Spearman \u03c1 (AUCell ~ cumul. AUROC)",
    y        = "Total AUCell change |ΔM3| + |ΔM27|"
  ) +
  theme_pub

save_fig(p_tf_cor, "fig3_tf_pk_correlation", width = 10, height = 7)

# -- F4: Rho distribution — Cmax vs cumulative AUROC comparison
p_rho_compare <- bind_rows(
  gene_pk_cor %>% transmute(item = gene, rho_AUCcum, rho_Cmax, type = "Gene"),
  tf_pk_cor   %>% transmute(item = tf,   rho_AUCcum, rho_Cmax, type = "TF")
) %>%
  pivot_longer(cols = c(rho_AUCcum, rho_Cmax),
               names_to = "pk_metric", values_to = "rho") %>%
  mutate(pk_metric = recode(pk_metric,
    rho_AUCcum = "Cumulative AUROC",
    rho_Cmax   = "Cmax"
  )) %>%
  ggplot(aes(x = rho, fill = pk_metric, colour = pk_metric)) +
  geom_density(alpha = 0.35, linewidth = 0.8) +
  geom_vline(xintercept = c(-0.9, -0.5, 0.5, 0.9),
             linetype = "dashed", colour = "grey50", linewidth = 0.5) +
  facet_wrap(~ type, ncol = 2) +
  scale_fill_manual(values   = c("Cumulative AUROC" = "#d95f02",
                                  "Cmax"             = "#2171b5"),
                    name = "PK metric") +
  scale_colour_manual(values = c("Cumulative AUROC" = "#d95f02",
                                  "Cmax"             = "#2171b5"),
                      name = "PK metric") +
  labs(
    title    = "Distribution of PK-expression Spearman correlations",
    subtitle = "Dashed lines: |ρ| = 0.5 (partial) and |ρ| = 0.9 (strong) thresholds",
    x        = "Spearman \u03c1",
    y        = "Density"
  ) +
  theme_pub

save_fig(p_rho_compare, "fig4_rho_distribution", width = 12, height = 6)

# -- F5: AutoML summary — R² and ΔR²
p_r2 <- summary_df %>%
  mutate(label = paste0(toupper(track), " \u2192", task)) %>%
  ggplot(aes(x = label, y = mean_R2, fill = label)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_errorbar(aes(ymin = mean_R2 - sd_R2,
                    ymax = mean_R2 + sd_R2), width = 0.25) +
  geom_point(aes(y = baseline_R2), shape = 4,
             colour = "black", size = 4, stroke = 1.5) +
  scale_fill_brewer(palette = "Set2") +
  labs(title    = "AutoML test R\u00b2 vs M0 baseline",
       subtitle = "\u00d7 = baseline-only R\u00b2  |  error bars = SD across noise replicates",
       x = NULL, y = expression(R^2~"(test)")) +
  theme_pub

p_delta <- summary_df %>%
  mutate(label = paste0(toupper(track), " \u2192 ", task)) %>%
  ggplot(aes(x = label, y = mean_delta_R2, fill = mean_delta_R2 > 0)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_errorbar(aes(ymin = mean_delta_R2 - sd_delta_R2,
                    ymax = mean_delta_R2 + sd_delta_R2), width = 0.25) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey40", linewidth = 0.7) +
  scale_fill_manual(values = c("TRUE" = "#2171b5", "FALSE" = "#cb181d")) +
  labs(title    = "PK contribution: \u0394R\u00b2 (AutoML \u2212 baseline)",
       subtitle = "Positive = drug exposure adds predictive value beyond M0",
       x = NULL, y = expression(Delta~R^2)) +
  theme_pub

p_automl <- p_r2 | p_delta
save_fig(p_automl, "fig5_automl_performance", width = 14, height = 6)

# -- F6: PK-sensitive vs PK-independent overlay on expression scatter
# (shows where PK-sensitive genes sit in the residual space)
gene_resid <- tibble(
  gene      = final_models[["gene_M3"]]$te_base$gene,
  resid_M3  = final_models[["gene_M3"]]$bl$test_resid,
  resid_M27 = final_models[["gene_M27"]]$bl$test_resid
) %>%
  left_join(gene_pk_cor %>% select(gene, rho_AUCcum, pk_class),
            by = "gene")

top_overlay <- gene_resid %>%
  filter(pk_class != "PK-independent") %>%
  mutate(score = abs(resid_M3) + abs(resid_M27)) %>%
  slice_max(score, n = 12, with_ties = FALSE)

p_overlay <- gene_resid %>%
  ggplot(aes(x = resid_M3, y = resid_M27, colour = pk_class)) +
  geom_point(alpha = 0.55, size = 2.0) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_text_repel(data = top_overlay, aes(label = gene),
                  size = BASE_FONT / 3.5,
                  max.overlaps = 20, show.legend = FALSE) +
  scale_colour_manual(values = pk_class_cols, name = "PK class") +
  labs(
    title    = "Gene residuals coloured by PK sensitivity class",
    subtitle = "PK-driven genes superimposed on M0-residualised expression space",
    x        = "M0-residualised expression at M3",
    y        = "M0-residualised expression at M27"
  ) +
  theme_pub

save_fig(p_overlay, "fig6_gene_residuals_pk_overlay", width = 10, height = 7)

# -- F7: Top PK-sensitive genes — expression trajectory lines
top_pk_genes <- gene_pk_cor %>%
  filter(pk_class %in% c("PK-driven (activated)", "PK-driven (suppressed)")) %>%
  slice_max(abs(rho_AUCcum), n = 12, with_ties = FALSE) %>%
  pull(gene)

if (length(top_pk_genes) > 0) {
  traj_df <- tibble(
    gene     = rep(top_pk_genes, each = 3),
    timepoint= rep(c("M0", "M3", "M27"), times = length(top_pk_genes)),
    expr     = c(rbind(expr_M0[top_pk_genes],
                       expr_M3[top_pk_genes],
                       expr_M27[top_pk_genes]))
  ) %>%
    left_join(gene_pk_cor %>% select(gene, rho_AUCcum, pk_class),
              by = "gene") %>%
    mutate(timepoint = factor(timepoint, levels = c("M0", "M3", "M27")))

  p_traj <- traj_df %>%
    ggplot(aes(x = timepoint, y = expr, group = gene, colour = pk_class)) +
    geom_line(linewidth = 1.0, alpha = 0.85) +
    geom_point(size = 3.0) +
    geom_text_repel(data = traj_df %>% filter(timepoint == "M27"),
                    aes(label = sprintf("%s (ρ=%.2f)", gene, rho_AUCcum)),
                    nudge_x = 0.15, hjust = 0,
                    size = 5, show.legend = FALSE) +
    scale_colour_manual(values = pk_class_cols, name = "PK class") +
    scale_x_discrete(expand = expansion(add = c(0.2, 0.8))) +
    labs(
      title    = "Expression trajectories of top PK-sensitive genes",
      #subtitle = "Genes with |ρ| \u2265 0.9 between expression and cumulative AUROC",
      x        = "Timepoint",
      y        = "Pseudo-bulk normalised expression"
    ) +
    theme_pub + 
    theme(
      axis.title       = element_text(size = 14),
      axis.text        = element_text(size = 12),
      legend.text      = element_text(size = 16),
      plot.title       = element_text(size = 14, face = "bold"),
      legend.title    = element_blank(),
      legend.position  = "bottom"
    )

  save_fig(p_traj, "fig7_pk_sensitive_gene_trajectories", width = 12, height = 9)
}

# -- F8: Combined main figure
combined <- (p_gene_cor | p_tf_cor) / (p_r2 | p_delta) +
  plot_annotation(
    title    = "PK-expression correlation analysis: ibrutinib resistance in CLL",
    subtitle = paste0(
      "Top: Spearman \u03c1 between expression/TF trajectories and cumulative AUROC. ",
      "Bottom: AutoML R\u00b2 performance."
    ),
    theme = theme(
      plot.title    = element_text(size = BASE_FONT + 3, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = BASE_FONT - 1, hjust = 0.5, colour = "grey40")
    )
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

save_fig(combined, "main_figure_combined", width = 16, height = 14)

## ============================================================
## 9. SUMMARY OUTPUT
## ============================================================

message("\n[8/8] Done.")
message("Outputs written to: ", OUTPUT_DIR)
message("")
message("Key files:")
message("  pk_summary_values.csv              — Cmax, AUC, cumulative AUROC at M3/M27")
message("  gene_pk_correlation.csv            — per-gene Spearman rho + PK class")
message("  tf_pk_correlation.csv              — per-TF Spearman rho + PK class")
message("  automl_summary.csv                 — AutoML R², ΔR², Spearman per track")
message("  fig1_pk_profile.svg/pdf            — PK profile with AUROC windows")
message("  fig2_gene_pk_correlation.svg/pdf   — gene rho scatter")
message("  fig3_tf_pk_correlation.svg/pdf     — TF rho scatter")
message("  fig4_rho_distribution.svg/pdf      — rho density: Cmax vs cumul. AUROC")
message("  fig5_automl_performance.svg/pdf    — R² and ΔR² comparison")
message("  fig6_gene_residuals_pk_overlay.svg — PK class overlay on residual scatter")
message("  fig7_pk_sensitive_gene_trajectories.svg — top PK-driven gene trajectories")
message("  main_figure_combined.svg/pdf       — combined 4-panel main figure")
