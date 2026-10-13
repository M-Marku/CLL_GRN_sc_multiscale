## ============================================================
## GRN-PK ibrutinib resistance study — full pipeline  (CORRECTED VERSION)
##
## Approach: PK-gene expression Spearman correlation
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
  library(svglite)
})

## ============================================================
## 0. CONFIGURATION
## ============================================================

SEURAT_PATH     <- "Data/cll_obj_tfAct_AUCell.rds"
PK_FILE         <- "Results/PK_results/ID6_AUC_comparaison_methodes_final.xlsx"
PK_SHEET        <- 1
DRIVER_M3_PATH  <- "Results/CellRank_results/driver_genes_cellfate3.csv"
DRIVER_M27_PATH <- "Results/CellRank_results/driver_genes_cellfate4.csv"
OUTPUT_DIR      <- "Results/ML_results/results_pk_correlation_revised"

# Clinical timepoints
M3_DAY          <- 90
M27_DAY         <- 822     # NOTE: the header comment says ~810 days; check which is right

# AutoML settings (used for the baseline R² validation block)
N_REPLICATES    <- 10
TEST_FRACTION   <- 0.20
AUTOML_TIME_SEC <- 60
H2O_MAX_MEM     <- "4G"
CV_FOLDS        <- 5
SEED            <- 42
NOISE_CV        <- 0.10    # jitter (CV) on the exposure value; see warning in section 7

# Trend classification (Spearman rho across THREE timepoints M0, M3, M27).
# With three distinct values rho can only be -1, -0.5, +0.5 or +1, so:
#   |rho| >= RHO_STRONG (0.9) means rho = +/-1, i.e. a monotonic trend over time;
#   |rho| <  RHO_WEAK   (0.5) can only occur with tied values.
# Cumulative AUC increases monotonically with time, so rho = +1 is identical to a
# correlation with time itself: exposure and time cannot be separated, and no
# p-value is possible with n = 3 (the smallest two-sided exact p is 1/3).
RHO_STRONG     <- 0.9
RHO_WEAK       <- 0.5
PK_LAB_UP      <- "Monotonic increase (rho = 1)"
PK_LAB_DOWN    <- "Monotonic decrease (rho = -1)"
PK_LAB_PARTIAL <- "Non-monotonic (|rho| = 0.5)"
PK_LAB_NONE    <- "No trend (|rho| < 0.5)"
# To keep the original names instead, set:
#   PK_LAB_UP <- "PK-driven (activated)";   PK_LAB_DOWN <- "PK-driven (suppressed)"
#   PK_LAB_PARTIAL <- "Partial PK correlation";   PK_LAB_NONE <- "PK-independent"

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
pk_all  <- pk_df   # use all rows (obs + imputed) for cumulative AUC

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
  # Cumulative AUC shading
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
  # AUC annotations
  annotate("text", x = M3_DAY / 2, y = y_max * 0.92,
           label = sprintf("Cumul. AUC\n0\u2192M3\n%.0f ng/mL\u00b7day", AUCcum_M3),
           size = BASE_FONT / 3.5, colour = "#2171b5",
           fontface = "bold", hjust = 0.5) +
  annotate("text", x = (M3_DAY + M27_DAY) / 2, y = y_max * 0.92,
           label = sprintf("Cumul. AUC\n0\u2192M27\n%.0f ng/mL\u00b7day", AUCcum_M27),
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
    title    = "Ibrutinib PK profile - observed vs simulated",
    subtitle = paste0(
      "Blue: observed daily median \u00b1 IQR (M1\u2013M6).  ",
      "Red: simulated daily mean Cmax (M0\u2013M27).\n",
      "Shaded: cumulative AUC windows (blue = 0\u2192M3; orange = 0\u2192M27)."
    )
  ) +
  theme_pub

save_fig(p_pk, "fig1_pk_profile", width = 15, height = 6)

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
## AUC trajectory (0, AUCcum_M3, AUCcum_M27).
## |rho| = 1: expression changes monotonically over M0 -> M3 -> M27, as
## cumulative AUC does. This cannot be told apart from a time trend, so it is
## a descriptive label, not evidence that exposure drives the change.
## ============================================================

message("\n[5/8] Computing PK-expression trajectory correlations...")

pk_traj <- c(0, AUCcum_M3, AUCcum_M27)   # cumulative AUC at M0, M3, M27
message(sprintf("  PK trajectory (cumul. AUC): M0=0  M3=%.0f  M27=%.0f ng/mL·day",
                AUCcum_M3, AUCcum_M27))

# Concentration (Cmax) is deliberately NOT tested against expression: trough
# concentration is flat from month 2, so Cmax at M3 and M27 is almost identical
# and rho would be decided by a negligible difference. Concentration is kept for
# the PK profile (fig1) and summary table only.

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
    # Spearman correlation: expression trajectory vs cumulative AUC
    rho_AUCcum = cor(
      c(expr_M0, expr_M3, expr_M27), pk_traj,
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
      abs(rho_AUCcum) >= RHO_STRONG & rho_AUCcum  > 0 ~ PK_LAB_UP,
      abs(rho_AUCcum) >= RHO_STRONG & rho_AUCcum  < 0 ~ PK_LAB_DOWN,
      abs(rho_AUCcum) <  RHO_WEAK                      ~ PK_LAB_NONE,
      TRUE                                              ~ PK_LAB_PARTIAL
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
    delta_M3  = auc_M3  - auc_M0,
    delta_M27 = auc_M27 - auc_M0,
    total_change = abs(delta_M3) + abs(delta_M27)
  ) %>%
  ungroup() %>%
  mutate(
    pk_class = case_when(
      abs(rho_AUCcum) >= RHO_STRONG & rho_AUCcum  > 0 ~ PK_LAB_UP,
      abs(rho_AUCcum) >= RHO_STRONG & rho_AUCcum  < 0 ~ PK_LAB_DOWN,
      abs(rho_AUCcum) <  RHO_WEAK                      ~ PK_LAB_NONE,
      TRUE                                              ~ PK_LAB_PARTIAL
    )
  ) %>%
  arrange(desc(abs(rho_AUCcum)))

write_csv(gene_pk_cor, file.path(OUTPUT_DIR, "gene_pk_correlation.csv"))
write_csv(tf_pk_cor,   file.path(OUTPUT_DIR, "tf_pk_correlation.csv"))

cat("\n--- Gene PK class distribution ---\n")
print(table(gene_pk_cor$pk_class))
cat("\n--- TF PK class distribution ---\n")
print(table(tf_pk_cor$pk_class))

# Diagnostic: with three timepoints only a few rho values exist
cat("\n--- Distinct rho values (genes) ---\n")
print(table(round(gene_pk_cor$rho_AUCcum, 2), useNA = "ifany"))
cat(sprintf("rho = +1: %.1f%% of genes | rho = -1: %.1f%% of genes\n",
            100 * mean(gene_pk_cor$rho_AUCcum >=  RHO_STRONG, na.rm = TRUE),
            100 * mean(gene_pk_cor$rho_AUCcum <= -RHO_STRONG, na.rm = TRUE)))
cat("Reference: for a random ordering of three values, 1/6 of genes give rho = +1 and 1/6 give rho = -1.\n")

cat("\n--- Top 15 genes by |rho| with cumulative AUC ---\n")
print(gene_pk_cor %>% select(gene, rho_AUCcum, pk_class) %>% head(15))

cat("\n--- Top 15 TFs by |rho| with cumulative AUC ---\n")
print(tf_pk_cor %>% select(tf, rho_AUCcum, pk_class) %>% head(15))

## ============================================================
## 7. AUTOML — BASELINE R² VALIDATION
##    Run AutoML with the real PK values from the Excel file.
##    This produces the ΔR² comparison figure for the manuscript.
## ============================================================

message("\n[6/8] h2o AutoML — nested comparison (corrected)...")

h2o.init(max_mem_size = H2O_MAX_MEM, nthreads = -1)
h2o.no_progress()

# Coefficient of determination on held-out genes: 1 - SSE/SST (can be < 0).
r2_score <- function(y, pred) {
  sst <- sum((y - mean(y))^2)
  if (!is.finite(sst) || sst == 0) return(NA_real_)
  1 - sum((y - pred)^2) / sst
}

# WARNING: exposure (AUCcum) is ONE value per timepoint, identical for
# every gene. Multiplying by (1 + N(0, NOISE_CV)) draws an independent random
# number per gene, i.e. a column unrelated to the gene's outcome. The "full"
# model therefore cannot beat the baseline-only model by design; dR2 ~ 0 is
# expected and is NOT evidence about PK effects.
# with_pk = FALSE gives the baseline-only feature set. `seed` makes the jitter
# differ between replicates (train and test use different seeds).
make_gene_X <- function(df, task, with_pk = TRUE, seed = SEED) {
  X <- df %>% transmute(expr_M0 = expr_M0, sd_M0 = sd_M0)
  if (with_pk) {
    pk_val <- if (task == "M3") AUCcum_M3 else AUCcum_M27
    n      <- nrow(df)
    set.seed(seed)
    X$AUCcum    <- pk_val * (1 + rnorm(n, 0, NOISE_CV))
    X$intx_cum  <- df$expr_M0 * X$AUCcum
  }
  X
}

make_tf_X <- function(df, task, with_pk = TRUE, seed = SEED) {
  X <- df %>% transmute(auc_M0 = auc_M0, auc_sd_M0 = sd_M0)
  if (with_pk) {
    pk_val <- if (task == "M3") AUCcum_M3 else AUCcum_M27
    n      <- nrow(df)
    set.seed(seed)
    X$AUCcum    <- pk_val * (1 + rnorm(n, 0, NOISE_CV))
    X$intx_cum  <- df$auc_M0 * X$AUCcum
  }
  X
}

# lm of M3/M27 on M0, fitted on TRAIN only; residuals define the AutoML target
remove_baseline <- function(train_df, test_df, target, predictor) {
  fit <- lm(as.formula(paste(target, "~", predictor)), data = train_df)
  list(
    train_resid = residuals(fit),
    test_resid  = test_df[[target]] - predict(fit, newdata = test_df),
    baseline_r2 = summary(fit)$r.squared,   # TRAIN R2 on RAW expression (reference only)
    fit         = fit
  )
}

# One AutoML fit; returns held-out predictions and metrics
fit_automl <- function(X_tr, y_tr, X_te, y_te, seed) {
  tr_h2o <- as.h2o(bind_cols(X_tr, resid = y_tr))
  te_h2o <- as.h2o(bind_cols(X_te, resid = y_te))
  aml <- h2o.automl(
    x = colnames(X_tr), y = "resid",
    training_frame     = tr_h2o,
    nfolds             = CV_FOLDS,
    max_runtime_secs   = AUTOML_TIME_SEC,
    seed               = seed,
    sort_metric        = "RMSE",
    stopping_metric    = "RMSE",
    stopping_rounds    = 3,
    keep_cross_validation_predictions = TRUE
  )
  pred <- as.vector(h2o.predict(aml@leader, te_h2o))
  out <- list(
    pred = pred,
    r2   = r2_score(y_te, pred),
    sp   = suppressWarnings(cor(y_te, pred, method = "spearman")),
    algo = str_extract(aml@leader@model_id,
                       "GBM|XGBoost|DRF|XRT|GLM|DeepLearning|StackedEnsemble")
  )
  h2o.rm(tr_h2o); h2o.rm(te_h2o)
  out
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
lm_ref       <- list()
final_models <- list()

for (rep in seq_len(N_REPLICATES)) {
  message(sprintf("\n  Replicate %d / %d", rep, N_REPLICATES))
  for (cfg in tasks_cfg) {
    label    <- paste0(cfg$track, "_", cfg$task_id)
    seed_rep <- SEED + rep * 100 +
                which(sapply(tasks_cfg, \(x)
                  x$task_id == cfg$task_id && x$track == cfg$track))

    tr <- cfg$train
    te <- cfg$test

    # Target = M0-residualised expression (lm fitted on train only)
    bl <- remove_baseline(tr, te, cfg$target, cfg$pred)
    tr$resid <- bl$train_resid
    te$resid <- bl$test_resid

    # Reference only: how much raw M3/M27 expression does M0 explain?
    if (rep == 1) {
      lm_ref[[label]] <- list(
        train = bl$baseline_r2,
        test  = r2_score(te[[cfg$target]], predict(bl$fit, newdata = te))
      )
    }

    # Nested comparison: SAME target, SAME genes, SAME seed
    res_base <- fit_automl(cfg$feat(tr, cfg$task_id, FALSE, seed_rep),       tr$resid,
                           cfg$feat(te, cfg$task_id, FALSE, seed_rep + 500), te$resid, seed_rep)
    res_full <- fit_automl(cfg$feat(tr, cfg$task_id, TRUE,  seed_rep),       tr$resid,
                           cfg$feat(te, cfg$task_id, TRUE,  seed_rep + 500), te$resid, seed_rep)

    all_metrics[[length(all_metrics) + 1]] <- tibble(
      track            = cfg$track, task = cfg$task_id,
      replicate        = rep, best_model = res_full$algo,
      R2_test          = res_full$r2,
      R2_baseline_only = res_base$r2,
      delta_R2         = res_full$r2 - res_base$r2,
      Spearman         = res_full$sp,
      lm_R2_raw_train  = lm_ref[[label]]$train,
      lm_R2_raw_test   = lm_ref[[label]]$test
    )

    if (rep == N_REPLICATES)
      final_models[[label]] <- list(pred = res_full$pred, y_te = te$resid,
                                    te_base = cfg$test, bl = bl)

    cat(sprintf("    [%s] rep%02d  %-16s  R2 base=%.3f  full=%.3f  dR2=%+.4f\n",
                label, rep,
                ifelse(is.na(res_full$algo), "unknown", res_full$algo),
                res_base$r2, res_full$r2, res_full$r2 - res_base$r2))
  }
}

metrics_df <- bind_rows(all_metrics)
summary_df <- metrics_df %>%
  group_by(track, task) %>%
  summarise(
    mean_R2         = mean(R2_test,          na.rm = TRUE),   # full model
    sd_R2           = sd(R2_test,            na.rm = TRUE),
    baseline_R2     = mean(R2_baseline_only, na.rm = TRUE),   # baseline-only AutoML, same target
    sd_baseline_R2  = sd(R2_baseline_only,   na.rm = TRUE),
    mean_delta_R2   = mean(delta_R2,         na.rm = TRUE),
    sd_delta_R2     = sd(delta_R2,           na.rm = TRUE),
    frac_delta_pos  = mean(delta_R2 > 0,     na.rm = TRUE),
    mean_Spearman   = mean(Spearman,         na.rm = TRUE),
    lm_R2_raw_train = first(lm_R2_raw_train),                 # M0 -> M3/M27 autocorrelation
    lm_R2_raw_test  = first(lm_R2_raw_test),
    .groups = "drop"
  )

write_csv(metrics_df,  file.path(OUTPUT_DIR, "automl_metrics.csv"))
write_csv(summary_df,  file.path(OUTPUT_DIR, "automl_summary.csv"))

h2o.shutdown(prompt = FALSE)

## ============================================================
## 8. FIGURES
## ============================================================

message("\n[7/8] Generating figures...")

# -- F2: PK-gene correlation scatter (cumulative AUC rho)
top_genes_pk <- gene_pk_cor %>%
  filter(pk_class != PK_LAB_NONE) %>%
  slice_max(abs(rho_AUCcum), n = 15, with_ties = FALSE)

pk_class_cols <- setNames(
  c("#d95f02", "#1f78b4", "#33a02c", "grey78"),
  c(PK_LAB_UP, PK_LAB_DOWN, PK_LAB_PARTIAL, PK_LAB_NONE)
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
    title    = "Gene expression trajectory vs cumulative AUC",
    #subtitle = "Spearman \u03c1 between expression (M0\u2192M3\u2192M27) and cumulative ibrutinib AUC",
    x        = "Spearman \u03c1 (expression ~ cumul. AUC)",
    y        = "Total expression change |ΔM3| + |ΔM27|"
  ) +
  theme_pub + 
  theme(
    axis.title       = element_text(size = 14),
    axis.text        = element_text(size = 13),
    legend.text      = element_text(size = 16),
    plot.title       = element_text(size = 14, face = "bold"),
    legend.title    = element_blank(),
    legend.position  = "bottom"
  )

save_fig(p_gene_cor, "fig2_gene_pk_correlation", width = 11, height = 7)

# -- F3: TF PK correlation scatter
top_tfs_pk <- tf_pk_cor %>%
  filter(pk_class != PK_LAB_NONE) %>%
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
    title    = "TF activity trajectory vs cumulative AUC",
    #subtitle = "Spearman \u03c1 between AUCell score (M0\u2192M3\u2192M27) and cumulative ibrutinib AUC",
    x        = "Spearman \u03c1 (AUCell ~ cumul. AUC)",
    y        = "Total AUCell change |ΔM3| + |ΔM27|"
  ) +
  theme_pub + 
  theme(
    axis.title       = element_text(size = 14),
    axis.text        = element_text(size = 13),
    legend.text      = element_text(size = 16),
    plot.title       = element_text(size = 14, face = "bold"),
    legend.title    = element_blank(),
    legend.position  = "bottom"
  )

save_fig(p_tf_cor, "fig3_tf_pk_correlation", width = 11, height = 7)

# (F4, the rho density comparison of Cmax vs cumulative AUC, was removed together with the Cmax test.)

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
       subtitle = "\u00d7 = baseline-only model R\u00b2 (same target, held-out genes)  |  error bars = SD across replicates",
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
  filter(pk_class != PK_LAB_NONE) %>%
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
    title    = "Gene residuals coloured by trend class",
    subtitle = "Trend genes (|rho| = 1) superimposed on M0-residualised expression space",
    x        = "M0-residualised expression at M3",
    y        = "M0-residualised expression at M27"
  ) +
  theme_pub

save_fig(p_overlay, "fig6_gene_residuals_pk_overlay", width = 10, height = 7)

# -- F6a/F6b: residual classes (M0-residualised change at M3 and M27), descriptive.
# NOTE: this classification used to live only in automl_grn_pk_h2o.R. It is now here.
# Unlike F6 above (held-out 20% of genes, coloured by trend class), these figures use ALL
# genes / regulons: the lm of M3 (or M27) on M0 is fitted on all features, and each feature is
# classified from its residuals (75th percentile of |residual| as cut-off, per timepoint).
classify_resid <- function(df, id, m0, m3, m27, q = 0.75) {
  r3  <- residuals(lm(reformulate(m0, m3),  data = df))
  r27 <- residuals(lm(reformulate(m0, m27), data = df))
  t3  <- quantile(abs(r3), q);  t27 <- quantile(abs(r27), q)
  tibble(id = df[[id]], resid_M3 = r3, resid_M27 = r27) %>%
    mutate(class = case_when(
      abs(resid_M27) > t27 & abs(resid_M3) <= t3                           ~ "M27-specific change",
      abs(resid_M3)  > t3  & abs(resid_M27) <= t27                         ~ "M3-specific change",
      abs(resid_M3)  > t3  & abs(resid_M27) > t27 &
        sign(resid_M3) == sign(resid_M27)                                  ~ "Consistent change",
      TRUE                                                                  ~ "Unchanged"))
}
resid_cols <- c("M27-specific change" = "#d95f02", "Consistent change" = "#1f78b4",
                "M3-specific change"  = "#33a02c", "Unchanged" = "grey78")

plot_resid <- function(res, title, xlab, ylab, n_label = 12, always_label = character()) {
  lab <- res %>% filter(class != "Unchanged") %>%
    slice_max(abs(resid_M3) + abs(resid_M27), n = n_label, with_ties = FALSE) %>%
    bind_rows(res %>% filter(id %in% always_label)) %>% distinct(id, .keep_all = TRUE)
  ggplot(res, aes(resid_M3, resid_M27, colour = class)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_text_repel(data = lab, aes(label = id), size = BASE_FONT / 3.5,
                    max.overlaps = 30, show.legend = FALSE) +
    scale_colour_manual(values = resid_cols, name = NULL) +
    labs(title    = title,
         subtitle = sprintf("n = %d; classes from M0-residuals (75th percentile cut-off)", nrow(res)),
         x = xlab, y = ylab) +
    theme_pub
}

gene_resid_cls <- classify_resid(base_gene_df, "gene", "expr_M0", "expr_M3", "expr_M27")
cat("\n--- Gene residual classes ---\n"); print(table(gene_resid_cls$class))
write_csv(gene_resid_cls, file.path(OUTPUT_DIR, "gene_residual_classes.csv"))
save_fig(plot_resid(gene_resid_cls,
                    "Genes: M0-residualised expression at M3 and M27",
                    "M0-residualised expression at M3", "M0-residualised expression at M27",
                    always_label = c("RGS2", "KLF2", "TCL1A", "IGLC2")),
         "fig6a_gene_residual_classes", width = 10, height = 7)

tf_resid_cls <- classify_resid(base_tf_df, "tf", "auc_M0", "auc_M3", "auc_M27")
cat("\n--- TF residual classes ---\n"); print(table(tf_resid_cls$class))
print(tf_resid_cls %>% filter(id %in% c("JUNB", "FOS", "JUN", "ATF4", "TP53")))   # check manuscript wording
write_csv(tf_resid_cls, file.path(OUTPUT_DIR, "tf_residual_classes.csv"))
save_fig(plot_resid(tf_resid_cls,
                    "TF regulons: M0-residualised activity at M3 and M27",
                    "M0-residualised AUCell at M3", "M0-residualised AUCell at M27",
                    always_label = c("JUNB", "FOS", "TP53")),
         "fig6b_tf_residual_classes", width = 10, height = 7)

# -- F7: Top PK-sensitive genes — expression trajectory lines
top_pk_genes <- gene_pk_cor %>%
  filter(pk_class %in% c(PK_LAB_UP, PK_LAB_DOWN)) %>%
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
      title    = "Expression trajectories of top genes by |rho| with cumulative AUC",
      #subtitle = "Genes with |ρ| \u2265 0.9 between expression and cumulative AUC",
      x        = "Timepoint",
      y        = "Pseudo-bulk normalised expression"
    ) +
    theme_pub + 
    theme(
      axis.title       = element_text(size = 14),
      axis.text        = element_text(size = 13),
      legend.text      = element_text(size = 16),
      plot.title       = element_text(size = 14, face = "bold"),
      legend.title    = element_blank(),
      legend.position  = "bottom"
    )

  save_fig(p_traj, "fig7_pk_sensitive_gene_trajectories", width = 10, height = 8.6)
}

# -- F8: Combined main figure
combined <- (p_gene_cor | p_tf_cor) / (p_r2 | p_delta) +
  plot_annotation(
    title    = "PK-expression correlation analysis: ibrutinib resistance in CLL",
    subtitle = paste0(
      "Top: Spearman \u03c1 between expression/TF trajectories and cumulative AUC. ",
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
message("  pk_summary_values.csv              — Cmax, AUC, cumulative AUC at M3/M27")
message("  gene_pk_correlation.csv            — per-gene Spearman rho + PK class")
message("  tf_pk_correlation.csv              — per-TF Spearman rho + PK class")
message("  automl_summary.csv                 — AutoML R², ΔR², Spearman per track")
message("  fig1_pk_profile.svg/pdf            — PK profile with AUC windows")
message("  fig2_gene_pk_correlation.svg/pdf   — gene rho scatter")
message("  fig3_tf_pk_correlation.svg/pdf     — TF rho scatter")
message("  fig5_automl_performance.svg/pdf    — R² and ΔR² comparison")
message("  fig6_gene_residuals_pk_overlay.svg — PK class overlay on residual scatter")
message("  fig7_pk_sensitive_gene_trajectories.svg — top genes by |rho| (trajectories)")
message("  main_figure_combined.svg/pdf       — combined 4-panel main figure")
