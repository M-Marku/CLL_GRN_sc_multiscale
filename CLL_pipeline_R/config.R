# =============================================================================
# config.R — shared paths, constants, and colour palettes for the CLL pipeline
# Source this file at the top of every numbered script: source("config.R")
# =============================================================================

# ── Directory structure (standardised — see README for the one path fix vs.
#    the original scripts: RSS_calculation.R used "../Plots_Results/" for one
#    figure; everything now goes under PLOTS_DIR like the rest) ---------------
DATA_DIR    <- "Data"
RESULTS_DIR <- "Results"
PLOTS_DIR   <- "Plots"

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS_DIR,   showWarnings = FALSE, recursive = TRUE)

# ── Key input file paths (edit here once, not per-script) --------------------
PBMC_FULL_RDATA   <- file.path(DATA_DIR, "PBMC/scRNAseq_obj_CLL.RData") 
CLL_ONLY_RDS      <- file.path(DATA_DIR, "CLL_only/cll_obj_tfAct_AUCell.rds")
AUCELL_CSV        <- file.path(RESULTS_DIR, "pyScenic_results/AUcell_CLL_only.csv")
FATE_PROBS_CSV    <- file.path(RESULTS_DIR, "CellRank_results/cell_matrix_fate_probabilities_light.csv")

# CopyKAT / phylogenetics — full PBMC object used as T-cell reference
PHYLO_OUT_DIR   <- "Results/Phylo_results"
CELLTYPE_COL    <- "cell_types"
REF_CELLTYPES   <- c("Tcells")
TIMEPOINT_COL   <- "time_point"

# ── Shared colour palette — timepoints (used in every script that plots by
#    timepoint: seurat UMAPs, AUCell heatmap/t-SNE, RSS, phylogenetics) -------
TIMEPOINT_COLORS <- c("M0" = "#66C2A5", "M3" = "#FC8D62", "M27" = "#8DA0CB")
TIMEPOINT_LEVELS <- c("M0", "M3", "M27")

# Clone colour palette (phylogenetics stage)
CLONE_COLORS <- c(
  "Clone_1"     = "#7F8C8D",
  "Clone_2"     = "#D4AC0D",
  "Clone_3"     = "#884EA0",
  "Diploid_ref" = "salmon"
)

set.seed(42)
