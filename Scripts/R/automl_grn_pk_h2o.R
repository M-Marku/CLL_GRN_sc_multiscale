## ============================================================
## AutoML pipeline — GRN-PK ibrutinib resistance study
## CLL cells: M0 (pre-treatment), M3 (response), M27 (relapse)
##
## AutoML engine: h2o::automl()
## Learners tried automatically:
##   - Gradient Boosting Machine (GBM)
##   - XGBoost
##   - Random Forest (DRF)
##   - Extremely Randomised Trees (XRT)
##   - ElasticNet (GLM)
##   - Deep Learning (MLP)
##   - Stacked Ensemble (best-of-family + all-models)
##
## Install:
##   install.packages("h2o")
##   install.packages(c("Seurat", "tidyverse", "ggrepel",
##                      "patchwork", "pheatmap", "viridis",
##                      "caret"))
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(ggrepel)
  library(patchwork)
  library(pheatmap)
  library(viridis)
  library(caret)    # for createDataPartition
  library(h2o)
})

## ============================================================
## 0. CONFIGURATION — edit paths and parameters before running
## ============================================================

SEURAT_PATH       <- "Data/CLL_only/cll_obj_tfAct_AUCell.rds"
DRIVER_M3_PATH    <- "Results/CLL_only/CellRank_results/driver_genes_cellfate3.csv"
DRIVER_M27_PATH   <- "Results/CLL_only/CellRank_results/driver_genes_cellfate4.csv"
AUCELL_ASSAY    <- "AUCell"
OUTPUT_DIR      <- "Results/CLL_only/ML_results/h2o_automl_conc"

#AUROC_M3        <- 42.7    # AUROC ibrutinib exposure at M3  (µg·h/mL)
#AUROC_M27       <- 38.1    # AUROC ibrutinib exposure at M27

NOISE_CV        <- 0.10    # Gaussian noise CV on AUROC (10%)
N_REPLICATES    <- 10      # noise replicates for stability
TEST_FRACTION   <- 0.20    # held-out test proportion
SEED            <- 42

# h2o AutoML budget (seconds per task per replicate).
# 60s is fine for exploration; use 300-600s for final results.
AUTOML_TIME_SEC <- 60

# h2o cluster memory (increase if you have RAM available)
H2O_MAX_MEM     <- "4G"

dir.create(OUTPUT_DIR, showWarnings = FALSE)


## ============================================================
## 1. INITIALISE h2o
## ============================================================

message("\n[Init] Starting h2o cluster...")
h2o.init(max_mem_size = H2O_MAX_MEM, nthreads = -1)  # -1 = all cores
h2o.no_progress()   # suppress verbose progress bars

## ============================================================
## 2. LOAD DATA
## ============================================================

message("\n[1/8] Loading data...")

seu <- readRDS(SEURAT_PATH)
DefaultAssay(seu) <- "RNA"

get_pseudobulk <- function(seurat_obj, tp) {
  cells <- colnames(seurat_obj)[seurat_obj$time_point == tp]
  rowMeans(GetAssayData(seurat_obj, layer = "data")[, cells, drop = FALSE])
}

get_pseudobulk_sd <- function(seurat_obj, tp) {
  cells <- colnames(seurat_obj)[seurat_obj$time_point == tp]
  apply(
    GetAssayData(seurat_obj, layer = "data")[, cells, drop = FALSE],
    1, sd, na.rm = TRUE
  )
}

expr_M0  <- get_pseudobulk(seu, "M0")
expr_M3  <- get_pseudobulk(seu, "M3")
expr_M27 <- get_pseudobulk(seu, "M27")
sd_M0    <- get_pseudobulk_sd(seu, "M0")

# AUCell TF matrix
DefaultAssay(seu) <- "AUCell"
auc_data <- GetAssayData(seu, layer = "data")  # TF × cell matrix

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

tfs_use <- rownames(auc_data)
message(sprintf("  TFs available: %d", length(tfs_use)))

# Restore default assay for gene expression
DefaultAssay(seu) <- "RNA"

## ============================================================
## 3. FEATURE MATRICES
## ============================================================

message("\n[2/8] Building feature matrices...")

drivers_M3  <- read_csv(DRIVER_M3_PATH) %>% as.data.frame()
drivers_M3 <- drivers_M3 %>% filter(drivers_M3$`3_corr`>0.1)
drivers_M3 <- drivers_M3[["...1"]]

drivers_M27 <- read_csv(DRIVER_M27_PATH) %>% as.data.frame()
drivers_M27 <- drivers_M27 %>% filter(drivers_M27$`4_corr`>0.1)
drivers_M27 <- drivers_M27[["...1"]]
genes_use   <- intersect(union(drivers_M3, drivers_M27), rownames(seu))

message(sprintf("  Driver genes: %d", length(genes_use)))

base_gene_df <- tibble(
  gene     = genes_use,
  expr_M0  = expr_M0[genes_use],
  expr_M3  = expr_M3[genes_use],
  expr_M27 = expr_M27[genes_use],
  sd_M0    = sd_M0[genes_use],
  lfc_M3   = expr_M3[genes_use] - expr_M0[genes_use],
  lfc_M27  = expr_M27[genes_use] - expr_M0[genes_use],
  in_M3    = genes_use %in% drivers_M3,
  in_M27   = genes_use %in% drivers_M27
) %>%
  filter(expr_M0 > 0) %>%
  mutate(expr_decile = ntile(expr_M0, 10))

base_tf_df <- tibble(
  tf          = tfs_use,
  auc_M0      = auc_M0[tfs_use],
  auc_M3      = auc_M3[tfs_use],
  auc_M27     = auc_M27[tfs_use],
  auc_sd_M0   = auc_sd_M0[tfs_use],
  auc_lfc_M3  = auc_M3[tfs_use] - auc_M0[tfs_use],
  auc_lfc_M27 = auc_M27[tfs_use] - auc_M0[tfs_use]
) %>%
  filter(auc_M0 > 0) %>%
  mutate(auc_decile = ntile(auc_M0, 10))

message(sprintf("  TFs retained: %d", nrow(base_tf_df)))

## ============================================================
## 4. TRAIN / TEST SPLIT   *** LEAKAGE FIREWALL ***
##    Fixed once before any modelling or baseline removal.
## ============================================================

message("\n[3/8] Train / test split (fixed, stratified by expression decile)...")

set.seed(SEED)
idx_gene_tr <- createDataPartition(
  base_gene_df$expr_decile, p = 1 - TEST_FRACTION, list = FALSE
)[, 1]
train_gene_base <- base_gene_df[ idx_gene_tr, ]
test_gene_base  <- base_gene_df[-idx_gene_tr, ]

set.seed(SEED + 1)
idx_tf_tr <- createDataPartition(
  base_tf_df$auc_decile, p = 1 - TEST_FRACTION, list = FALSE
)[, 1]
train_tf_base <- base_tf_df[ idx_tf_tr, ]
test_tf_base  <- base_tf_df[-idx_tf_tr, ]

message(sprintf("  Gene — train: %d  test: %d",
                nrow(train_gene_base), nrow(test_gene_base)))
message(sprintf("  TF   — train: %d  test: %d",
                nrow(train_tf_base),  nrow(test_tf_base)))

## ============================================================
## 5. HELPER FUNCTIONS
## ============================================================

# Gaussian noise on AUROC values (addresses scalar PK degeneracy)
add_noise <- function(df, noise_cv, seed) {
  set.seed(seed)
  n <- nrow(df)
  df$AUROC_M3  <- AUROC_M3  * (1 + rnorm(n, 0, noise_cv))
  df$AUROC_M27 <- AUROC_M27 * (1 + rnorm(n, 0, noise_cv))
  df
}

# Baseline autocorrelation removal
# lm fitted on TRAIN only, residuals computed for both splits
remove_baseline <- function(train_df, test_df, target, predictor) {
  fm   <- as.formula(paste(target, "~", predictor))
  fit  <- lm(fm, data = train_df)
  list(
    train_resid = residuals(fit),
    test_resid  = test_df[[target]] - predict(fit, newdata = test_df),
    baseline_r2 = summary(fit)$r.squared,
    fit         = fit
  )
}

# Feature matrix for gene track
# Feature matrix for gene track — AUROC is task-specific
make_gene_X <- function(df, task) {
  auroc_col <- if (task == "M3") df$AUROC_M3 else df$AUROC_M27
  df %>% transmute(
    expr_M0  = expr_M0,
    sd_M0    = sd_M0,
    AUROC    = auroc_col,          # only the relevant timepoint AUROC
    intx     = expr_M0 * auroc_col # M0 × drug interaction
  )
}

# Feature matrix for TF track — same logic
make_tf_X <- function(df, task) {
  auroc_col <- if (task == "M3") df$AUROC_M3 else df$AUROC_M27
  df %>% transmute(
    auc_M0    = auc_M0,
    auc_sd_M0 = auc_sd_M0,
    AUROC     = auroc_col,
    intx      = auc_M0 * auroc_col
  )
}

# Run h2o AutoML and return key results
run_h2o_automl <- function(train_X, train_y, test_X, test_y,
                            task_id, rep, seed) {

  # Assemble h2o data frames
  tr_h2o <- as.h2o(bind_cols(train_X, resid = train_y))
  te_h2o <- as.h2o(bind_cols(test_X,  resid = test_y))

  predictors <- colnames(train_X)
  response   <- "resid"

  # Run AutoML
  aml <- h2o.automl(
    x                  = predictors,
    y                  = response,
    training_frame     = tr_h2o,
    nfolds             = 5,              # inner CV folds
    max_runtime_secs   = AUTOML_TIME_SEC,
    seed               = seed,
    sort_metric        = "RMSE",         # internal ranking metric
    stopping_metric    = "RMSE",
    stopping_rounds    = 3,
    keep_cross_validation_predictions = TRUE,   # required for stacking
    exclude_algos      = NULL            # try everything available
  )

  # Extract leaderboard
  lb <- as.data.frame(h2o.get_leaderboard(aml, extra_columns = "ALL")) %>%
    select(model_id, rmse, mse, mae, rmsle) %>%
    mutate(task = task_id, replicate = rep)

  # Best model predictions on test set
  best_model <- aml@leader
  pred_te    <- as.vector(h2o.predict(best_model, te_h2o))

  # Test metrics
  r2   <- if (sd(test_y) > 0) cor(test_y, pred_te)^2 else NA_real_
  sp   <- cor(test_y, pred_te, method = "spearman")
  rmse <- sqrt(mean((test_y - pred_te)^2))

  # Variable importance from best model (not all model types support it)
  varimp <- tryCatch(
    as.data.frame(h2o.varimp(best_model)) %>%
      mutate(task = task_id, replicate = rep, model = best_model@model_id),
    error = function(e) NULL
  )

  # Clean up h2o frames to free memory
  h2o.rm(tr_h2o)
  h2o.rm(te_h2o)

  list(
    leaderboard  = lb,
    best_model_id= best_model@model_id,
    pred         = pred_te,
    r2           = r2,
    spearman     = sp,
    rmse         = rmse,
    varimp       = varimp
  )
}

## ============================================================
## 6. TASK CONFIGURATION
## ============================================================

tasks_cfg <- list(
  list(track = "gene", task_id = "M3",
       target = "expr_M3",  predictor = "expr_M0",
       train_base = train_gene_base, test_base = test_gene_base,
       feat_fn = make_gene_X),
  list(track = "gene", task_id = "M27",
       target = "expr_M27", predictor = "expr_M0",
       train_base = train_gene_base, test_base = test_gene_base,
       feat_fn = make_gene_X),
  list(track = "tf",   task_id = "M3",
       target = "auc_M3",   predictor = "auc_M0",
       train_base = train_tf_base, test_base = test_tf_base,
       feat_fn = make_tf_X),
  list(track = "tf",   task_id = "M27",
       target = "auc_M27",  predictor = "auc_M0",
       train_base = train_tf_base, test_base = test_tf_base,
       feat_fn = make_tf_X)
)

## ============================================================
## 7. MAIN LOOP: replicates × tasks
## ============================================================

message(sprintf(
  "\n[4/8] h2o AutoML — %d replicates × %d tasks (%d s budget each)...",
  N_REPLICATES, length(tasks_cfg), AUTOML_TIME_SEC
))

all_metrics   <- list()
all_lbs       <- list()
all_varimp    <- list()
baseline_r2   <- list()   # stored once from replicate 1
final_results <- list()   # predictions from last replicate for classification

for (rep in seq_len(N_REPLICATES)) {
  message(sprintf("\n  === Replicate %d / %d ===", rep, N_REPLICATES))
  
  for (cfg in tasks_cfg) {
    label    <- paste0(cfg$track, "_", cfg$task_id)
    seed_rep <- SEED + rep * 100 + which(sapply(tasks_cfg, \(x)
                                                x$task_id == cfg$task_id && x$track == cfg$track))
    
    # Add AUROC noise independently to train and test
    tr <- add_noise(cfg$train_base, NOISE_CV, seed_rep)
    te <- add_noise(cfg$test_base,  NOISE_CV, seed_rep + 500)
    
    # Baseline removal (fitted on train only, applied to test)
    bl <- remove_baseline(tr, te, cfg$target, cfg$predictor)
    tr$resid <- bl$train_resid
    te$resid <- bl$test_resid
    
    if (rep == 1) baseline_r2[[label]] <- bl$baseline_r2
    
    # Feature matrices — AUROC is task-specific (no cross-timepoint leakage)
    X_tr <- cfg$feat_fn(tr, cfg$task_id)
    X_te <- cfg$feat_fn(te, cfg$task_id)
    
    # Run h2o AutoML
    res <- run_h2o_automl(
      train_X = X_tr, train_y = tr$resid,
      test_X  = X_te, test_y  = te$resid,
      task_id = label, rep = rep, seed = seed_rep
    )
    
    # Collect metrics
    all_metrics[[length(all_metrics) + 1]] <- tibble(
      track        = cfg$track,
      task         = cfg$task_id,
      replicate    = rep,
      best_model   = res$best_model_id,
      R2_test      = res$r2,
      Spearman     = res$spearman,
      RMSE_test    = res$rmse,
      baseline_R2  = baseline_r2[[label]],
      delta_R2     = res$r2 - baseline_r2[[label]]
    )
    
    if (!is.null(res$leaderboard))
      all_lbs[[length(all_lbs) + 1]] <- res$leaderboard
    
    if (!is.null(res$varimp))
      all_varimp[[length(all_varimp) + 1]] <- res$varimp
    
    # Keep last replicate results for downstream classification
    if (rep == N_REPLICATES) {
      final_results[[label]] <- list(
        pred    = res$pred,
        y_te    = te$resid,
        te_base = cfg$test_base,
        bl      = bl,
        X_te    = X_te
      )
    }
    
    # Extract algorithm family from model id for display
    algo <- str_extract(res$best_model_id,
                        "GBM|XGBoost|DRF|XRT|GLM|DeepLearning|StackedEnsemble")
    cat(sprintf("    [%s] rep%02d  best: %-16s  R²=%.3f  ΔR²=%+.4f  ρ=%.3f\n",
                label, rep,
                ifelse(is.na(algo), "unknown", algo),
                res$r2,
                res$r2 - baseline_r2[[label]],
                res$spearman))
  }
}

## ============================================================
## 8. AGGREGATE & SAVE RESULTS
## ============================================================

message("\n[5/8] Aggregating results...")

metrics_df <- bind_rows(all_metrics)

summary_df <- metrics_df %>%
  group_by(track, task) %>%
  summarise(
    mean_R2        = mean(R2_test,    na.rm = TRUE),
    sd_R2          = sd(R2_test,      na.rm = TRUE),
    mean_delta_R2  = mean(delta_R2,   na.rm = TRUE),
    sd_delta_R2    = sd(delta_R2,     na.rm = TRUE),
    mean_Spearman  = mean(Spearman,   na.rm = TRUE),
    baseline_R2    = first(baseline_R2),
    .groups = "drop"
  )

# Model family frequency (strip h2o model ID suffix to get algorithm name)
model_freq <- metrics_df %>%
  mutate(algo = str_extract(best_model,
    "GBM|XGBoost|DRF|XRT|GLM|DeepLearning|StackedEnsemble")) %>%
  mutate(algo = replace_na(algo, "Other")) %>%
  count(track, task, algo) %>%
  group_by(track, task) %>%
  mutate(pct = n / sum(n) * 100) %>%
  arrange(track, task, desc(n))

# Leaderboard across all tasks and replicates
lb_all <- bind_rows(all_lbs) %>%
  mutate(algo = str_extract(model_id,
    "GBM|XGBoost|DRF|XRT|GLM|DeepLearning|StackedEnsemble"))

# Variable importance averaged across replicates
varimp_avg <- bind_rows(all_varimp) %>%
  group_by(task, variable) %>%
  summarise(mean_importance = mean(scaled_importance, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(task, desc(mean_importance))

write_csv(metrics_df,   file.path(OUTPUT_DIR, "automl_metrics_all.csv"))
write_csv(summary_df,   file.path(OUTPUT_DIR, "automl_summary.csv"))
write_csv(model_freq,   file.path(OUTPUT_DIR, "automl_model_freq.csv"))
write_csv(lb_all,       file.path(OUTPUT_DIR, "automl_leaderboard_all.csv"))
write_csv(varimp_avg,   file.path(OUTPUT_DIR, "automl_varimp_avg.csv"))

cat("\n===== SUMMARY =====\n")
print(summary_df)
cat("\n===== MODEL SELECTION FREQUENCY =====\n")
print(model_freq)

## ============================================================
## 9. GENE & TF CLASSIFICATION
## ============================================================

message("\n[6/8] Gene and TF classification...")

classify_items <- function(resid_M3, resid_M27) {
  q75_M3  <- quantile(abs(resid_M3),  0.75, na.rm = TRUE)
  q75_M27 <- quantile(abs(resid_M27), 0.75, na.rm = TRUE)
  case_when(
    abs(resid_M27) > q75_M27 & abs(resid_M3) <= q75_M3  ~ "Acquired resistance",
    abs(resid_M3)  > q75_M3  & abs(resid_M27) > q75_M27 &
      sign(resid_M3) == sign(resid_M27)                  ~ "Intrinsic",
    abs(resid_M3)  > q75_M3  & abs(resid_M27) <= q75_M27 ~ "Response-specific",
    TRUE                                                   ~ "Neutral"
  )
}

gene_class_df <- final_results[["gene_M3"]]$te_base %>%
  select(gene, expr_M0, expr_M3, expr_M27) %>%
  mutate(
    resid_M3  = final_results[["gene_M3"]]$bl$test_resid,
    resid_M27 = final_results[["gene_M27"]]$bl$test_resid,
    class     = classify_items(resid_M3, resid_M27)
  )

tf_class_df <- final_results[["tf_M3"]]$te_base %>%
  select(tf, auc_M0, auc_M3, auc_M27) %>%
  mutate(
    resid_M3  = final_results[["tf_M3"]]$bl$test_resid,
    resid_M27 = final_results[["tf_M27"]]$bl$test_resid,
    class     = classify_items(resid_M3, resid_M27)
  )

write_csv(gene_class_df, file.path(OUTPUT_DIR, "gene_classification.csv"))
write_csv(tf_class_df,   file.path(OUTPUT_DIR, "tf_classification.csv"))

cat("\n--- Gene class distribution ---\n"); print(table(gene_class_df$class))
cat("\n--- TF class distribution ---\n");   print(table(tf_class_df$class))

## ============================================================
## 10. FIGURES
## ============================================================

message("\n[7/8] Generating figures...")

theme_pub <- theme_minimal(base_size = 11) +
  theme(strip.text      = element_text(face = "bold"),
        legend.position = "bottom",
        plot.title      = element_text(face = "bold"))

class_cols <- c(
  "Acquired resistance" = "#d95f02",
  "Intrinsic"          = "#1f78b4",
  "Response-specific"  = "#33a02c",
  "Neutral"            = "grey78"
)

# -- F1: Model selection frequency heatmap
p_freq <- model_freq %>%
  mutate(label = paste0(toupper(track), " → M", task)) %>%
  ggplot(aes(x = label, y = algo, fill = pct)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.0f%%", pct)), size = 3.2) +
  scale_fill_viridis_c(name = "% selected", limits = c(0, 100),
                       option = "plasma") +
  labs(title = "AutoML — winning algorithm per task",
       x = NULL, y = NULL) +
  theme_pub

ggsave(file.path(OUTPUT_DIR, "fig1_model_freq_heatmap.svg"),
       p_freq, width = 8, height = 4.5, device = svg)

# -- F2: Test R² vs baseline
p_r2 <- summary_df %>%
  mutate(label = paste0(toupper(track), " → ", task)) %>%
  ggplot(aes(x = label, y = mean_R2, fill = label)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_errorbar(aes(ymin = mean_R2 - sd_R2,
                    ymax = mean_R2 + sd_R2), width = 0.2) +
  geom_point(aes(y = baseline_R2), shape = 4,
             colour = "black", size = 3.5, stroke = 1.3) +
  scale_fill_brewer(palette = "Set2") +
  labs(title    = "AutoML test R² vs baseline (× = M0 autocorrelation R²)",
       subtitle = "Error bars = SD across noise replicates",
       x = NULL, y = expression(R^2~"(test)")) +
  theme_pub + 
  theme(
    axis.title       = element_text(size = 14),
    axis.text        = element_text(size = 14),
    legend.text      = element_text(size = 16),
    plot.title       = element_text(size = 16, face = "bold"),
    legend.title    = element_blank(),
    legend.position  = "bottom"
  )

ggsave(file.path(OUTPUT_DIR, "fig2_r2_comparison.svg"),
       p_r2, width = 7, height = 5, device = svg)

# -- F3: AUROC ΔR²
p_delta <- summary_df %>%
  mutate(label = paste0(toupper(track), " → ", task)) %>%
  ggplot(aes(x = label, y = mean_delta_R2, fill = mean_delta_R2 > 0)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_errorbar(aes(ymin = mean_delta_R2 - sd_delta_R2,
                    ymax = mean_delta_R2 + sd_delta_R2), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = c("TRUE" = "#2171b5", "FALSE" = "#cb181d")) +
  labs(title    = "AUROC drug exposure contribution (ΔR²)",
       subtitle = "AutoML R² minus baseline-only R²",
       x = NULL, y = expression(Delta~R^2)) +
  theme_pub + 
  theme(
    axis.title       = element_text(size = 14),
    axis.text        = element_text(size = 14),
    legend.text      = element_text(size = 16),
    plot.title       = element_text(size = 16, face = "bold"),
    legend.title    = element_blank(),
    legend.position  = "bottom"
  )

ggsave(file.path(OUTPUT_DIR, "fig3_delta_r2.svg"),
       p_delta, width = 8, height = 5, device = svg)

# -- F4: Gene classification scatter
top_genes <- gene_class_df %>%
  filter(class != "Neutral") %>%
  mutate(score = abs(resid_M3) + abs(resid_M27)) %>%
  slice_max(score, n = 14, with_ties = FALSE)

p_gene <- gene_class_df %>%
  ggplot(aes(x = resid_M3, y = resid_M27, colour = class)) +
  geom_point(alpha = 0.55, size = 1.6) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_text_repel(data = top_genes, aes(label = gene),
                  size = 5, max.overlaps = 20, show.legend = FALSE) +
  scale_colour_manual(values = class_cols, name = NULL) +
  labs(title = "Gene mechanistic classification",
       x = "M0-residualised expression at M3",
       y = "M0-residualised expression at M27") +
  theme_pub + 
  theme(
    axis.title       = element_text(size = 14),
    axis.text        = element_text(size = 14),
    legend.text      = element_text(size = 16),
    plot.title       = element_text(size = 16, face = "bold"),
    legend.title    = element_blank(),
    legend.position  = "bottom"
  )

ggsave(file.path(OUTPUT_DIR, "fig4_gene_classification.svg"),
       p_gene, width = 8, height = 7, device = svg)

# -- F5: TF classification scatter
top_tfs <- tf_class_df %>%
  filter(class != "Neutral") %>%
  mutate(score = abs(resid_M3) + abs(resid_M27)) %>%
  slice_max(score, n = 14, with_ties = FALSE)

p_tf <- tf_class_df %>%
  ggplot(aes(x = resid_M3, y = resid_M27, colour = class)) +
  geom_point(alpha = 0.55, size = 1.6) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_text_repel(data = top_tfs, aes(label = tf),
                  size = 5, max.overlaps = 20, show.legend = FALSE) +
  scale_colour_manual(values = class_cols, name = NULL) +
  labs(title = "TF mechanistic classification (AUCell residuals)",
       x = "AUCell residual at M3",
       y = "AUCell residual at M27") +
  theme_pub + 
  theme(
    axis.title       = element_text(size = 14),
    axis.text        = element_text(size = 14),
    legend.text      = element_text(size = 16),
    plot.title       = element_text(size = 16, face = "bold"),
    legend.title    = element_blank(),
    legend.position  = "bottom"
  )

ggsave(file.path(OUTPUT_DIR, "fig5_tf_classification.svg"),
       p_tf, width = 8, height = 7, device = svg)

# -- F6: Variable importance (averaged across replicates, gene M27 track)
p_varimp <- varimp_avg %>%
  filter(str_detect(task, "gene_M27")) %>%
  slice_max(mean_importance, n = 10) %>%
  ggplot(aes(x = reorder(variable, mean_importance),
             y = mean_importance)) +
  geom_col(fill = "#2171b5", width = 0.6) +
  coord_flip() +
  labs(title    = "Variable importance - gene M27 track",
       subtitle = "Mean scaled importance across noise replicates",
       x = NULL, y = "Scaled importance") +
  theme_pub + 
  theme(
    axis.title       = element_text(size = 14),
    axis.text        = element_text(size = 14),
    legend.text      = element_text(size = 16),
    plot.title       = element_text(size = 16, face = "bold"),
    legend.title    = element_blank(),
    legend.position  = "bottom"
  )

ggsave(file.path(OUTPUT_DIR, "fig6_varimp_geneM27.svg"),
       p_varimp, width = 6, height = 4.5, device = svg)

# -- F7: Leaderboard RMSE distribution (box per algorithm family)
p_lb <- lb_all %>%
  filter(!is.na(algo)) %>%
  ggplot(aes(x = reorder(algo, rmse, median), y = rmse, fill = algo)) +
  geom_boxplot(outlier.size = 0.8, show.legend = FALSE) +
  facet_wrap(~ task, scales = "free_y", ncol = 2) +
  coord_flip() +
  scale_fill_brewer(palette = "Set3") +
  labs(title    = "AutoML leaderboard RMSE by algorithm family",
       subtitle = "Across all tasks and noise replicates",
       x = NULL, y = "CV RMSE") +
  theme_pub

ggsave(file.path(OUTPUT_DIR, "fig7_leaderboard_rmse.svg"),
       p_lb, width = 10, height = 7, device = svg)

# -- F8: Pred vs observed (gene M27, last replicate)
fr_g27 <- final_results[["gene_M27"]]
p_po <- tibble(obs = fr_g27$y_te, pred = fr_g27$pred) %>%
  ggplot(aes(obs, pred)) +
  geom_point(alpha = 0.45, size = 1.4, colour = "#2171b5") +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey40") +
  geom_smooth(method = "lm", se = FALSE,
              colour = "#d95f02", linewidth = 0.8) +
  labs(title    = "Predicted vs observed — gene M27 (best model)",
       x = "Observed M0-residualised expression",
       y = "Predicted M0-residualised expression") +
  theme_pub

ggsave(file.path(OUTPUT_DIR, "fig8_pred_vs_obs.svg"),
       p_po, width = 6, height = 5, device = svg)

## ============================================================
## 11. COMBINED MAIN FIGURE
## ============================================================

message("\n[8/8] Saving combined main figure...")

combined <- (p_r2 | p_delta) / (p_gene | p_tf) +
  plot_annotation(
    title    = "h2o AutoML: PK–RNAseq integration in CLL ibrutinib resistance",
    subtitle = "Top: predictive performance. Bottom: gene and TF mechanistic classification.",
    theme    = theme(
      plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, colour = "grey40")
    )
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(OUTPUT_DIR, "main_figure_combined.svg"),
       combined, width = 14, height = 12, device = svg)

## ============================================================
## 12. SHUT DOWN h2o
## ============================================================

h2o.shutdown(prompt = FALSE)

message("\n===== DONE =====")
message("Outputs written to: ", OUTPUT_DIR)
message("")
message("Key files:")
message("  automl_summary.csv          — mean R², ΔR², Spearman per task")
message("  automl_model_freq.csv       — winning algorithm frequency")
message("  automl_leaderboard_all.csv  — full leaderboard, all tasks/replicates")
message("  automl_varimp_avg.csv       — variable importance averaged across replicates")
message("  gene_classification.csv     — gene mechanistic classes")
message("  tf_classification.csv       — TF mechanistic classes")
message("  main_figure_combined.pdf    — main text figure (R² + ΔR² + scatters)")
message("  fig1_model_freq_heatmap.pdf — which algorithm won per task")
message("  fig7_leaderboard_rmse.pdf   — RMSE distribution by algorithm family")

