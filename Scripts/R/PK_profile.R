## ============================================================
## PK profile plot + cumulative AUROC calculation
## Ibrutinib PK data: hourly measurements, observed until M6,
## simulated until M27.
##
## Cumulative AUROC definition (correct for resistance question):
##   AUROC_M3  = integral of Cp_IBRU_pred from day 0 to day 90
##   AUROC_M27 = integral of Cp_IBRU_pred from day 0 to day 810
## M27 > M3 by construction (cumulative exposure).
##
## Columns used from Excel:
##   time_abs_d    : absolute time in days
##   Cp_IBRU_pred  : predicted/simulated ibrutinib concentration (ng/mL)
##   Cp_IBRU_obs   : observed ibrutinib concentration (ng/mL)
##   imputed       : flag for imputed rows (0 = real, 1 = imputed)
##
## Output:
##   pk_profile.pdf / .png  — publication-ready PK figure
##   Console: AUROC_M3 and AUROC_M27 values to plug into the ML pipeline
## ============================================================

library(readxl)
library(ggplot2)
library(dplyr)

## ============================================================
## 0. CONFIGURATION
## ============================================================

PK_FILE    <- "Results/CLL_only/PK_results/ID6_AUC_comparaison_methodes_final.xlsx"
PK_SHEET   <- 1
OUTPUT_DIR <- "Results/CLL_only/PK_results/"   # same output dir as the main pipeline

# Clinical timepoints in days
M3_DAY  <- 90    # 3 months
M6_DAY  <- 180   # 6 months — last observed concentration
M27_DAY <- 810   # 27 months — last simulated timepoint

dir.create(OUTPUT_DIR, showWarnings = FALSE)

## ============================================================
## 1. LOAD AND CLEAN
## ============================================================

message("\n[1/4] Loading PK data...")

pk_raw <- read_excel(PK_FILE, sheet = PK_SHEET)

message(sprintf("  Rows loaded: %d", nrow(pk_raw)))
message(sprintf("  Columns: %s", paste(colnames(pk_raw), collapse = ", ")))

pk_df <- pk_raw %>%
  rename(
    time_d = day,
    pred   = Cmax_IBRU,
    obs    = Ctrough_IBRU
  ) %>%
  filter(time_d >= 0, time_d <= M27_DAY) %>%
  #filter(imputed == FALSE) %>%
  filter(!is.na(pred) | !is.na(obs))

message(sprintf("  Rows after cleaning: %d", nrow(pk_df)))
message(sprintf("  Time range: %.1f – %.1f days", min(pk_df$time_d), max(pk_df$time_d)))

## ============================================================
## 2. CUMULATIVE AUROC (trapezoidal integration on predicted curve)
##
##   Cumulative from day 0 in both cases — this is the correct
##   definition for the question "does accumulated drug exposure
##   drive resistance?". AUROC_M27 > AUROC_M3 by construction.
## ============================================================

message("\n[2/4] Computing cumulative AUROC...")

pk_pred_full <- pk_df %>%
  filter(!is.na(pred)) %>%
  arrange(time_d)

# Cumulative AUROC 0 → M3 (day 0 to day 90)
auroc_M3 <- pk_pred_full %>%
  filter(time_d <= M3_DAY) %>%
  summarise(
    auroc = sum(diff(time_d) *
                  (head(pred, -1) + tail(pred, -1)) / 2)
  ) %>%
  pull(auroc)

# Cumulative AUROC 0 → M27 (day 0 to day 810)
auroc_M27 <- pk_pred_full %>%
  filter(time_d <= M27_DAY) %>%
  summarise(
    auroc = sum(diff(time_d) *
                  (head(pred, -1) + tail(pred, -1)) / 2)
  ) %>%
  pull(auroc)

# Sanity check: row counts per window
n_M3  <- pk_pred_full %>% filter(time_d <= M3_DAY)  %>% nrow()
n_M27 <- pk_pred_full %>% filter(time_d <= M27_DAY) %>% nrow()

message(sprintf("  Rows used for AUROC_M3  (0–%.0fd):  %d", M3_DAY,  n_M3))
message(sprintf("  Rows used for AUROC_M27 (0–%.0fd): %d", M27_DAY, n_M27))
message(sprintf("  Cumulative AUROC M3  (0→%.0fd):  %.1f ng/mL·day", M3_DAY,  auroc_M3))
message(sprintf("  Cumulative AUROC M27 (0→%.0fd): %.1f ng/mL·day", M27_DAY, auroc_M27))
message(sprintf("  Ratio M27/M3: %.2f (expected ~%.1f given window ratio)",
                auroc_M27 / auroc_M3, M27_DAY / M3_DAY))

message("\n  *** Plug these into the ML pipeline: ***")
message(sprintf("  AUROC_M3  <- %.1f", auroc_M3))
message(sprintf("  AUROC_M27 <- %.1f", auroc_M27))

## ============================================================
## 3. PREPARE PLOT DATA
## ============================================================

message("\n[3/4] Preparing plot data...")

# Observed: summarise to daily median + IQR (M0 to M6 only)
pk_daily_obs <- pk_df %>%
  filter(!is.na(obs), obs >= 0, time_d <= M6_DAY) %>%
  group_by(day = floor(time_d)) %>%
  summarise(
    median_obs = median(obs,                       na.rm = TRUE),
    q25_obs    = quantile(obs, 0.25,               na.rm = TRUE),
    q75_obs    = quantile(obs, 0.75,               na.rm = TRUE),
    .groups = "drop"
  )

# Simulated: daily trough (min per day) for clean trend line
pk_pred_daily <- pk_pred_full %>%
  mutate(day_round = round(time_d)) %>%
  group_by(day_round) %>%
  summarise(pred = mean(pred, na.rm = TRUE), .groups = "drop")

# Y axis ceiling: 99th percentile of observed + 20% headroom
y_max <- quantile(pk_df$obs, 0.99, na.rm = TRUE) * 1.20

message(sprintf("  Observed days available: %d (M0–M6)",  nrow(pk_daily_obs)))
message(sprintf("  Simulated days available: %d (M0–M27)", nrow(pk_pred_daily)))
message(sprintf("  Y axis ceiling: %.1f ng/mL", y_max))

## ============================================================
## 4. PLOT
## ============================================================

message("\n[4/4] Generating figure...")

p_pk <- ggplot() +
  
  ## -- Cumulative AUROC shading
  # M27 window (full, behind M3 shading)
  annotate("rect",
           xmin = 0, xmax = M27_DAY, ymin = 0, ymax = Inf,
           fill = "#d95f02", alpha = 0.05) +
  # M3 window (on top, darker blue)
  annotate("rect",
           xmin = 0, xmax = M3_DAY,  ymin = 0, ymax = Inf,
           fill = "#2171b5", alpha = 0.10) +
  
  ## -- Observed IQR ribbon (M0–M6)
  geom_ribbon(
    data = pk_daily_obs,
    aes(x = day, ymin = q25_obs, ymax = q75_obs),
    fill = "#2171b5", alpha = 0.22
  ) +
  
  ## -- Observed daily median points (M0–M6)
  geom_point(
    data = pk_daily_obs,
    aes(x = day, y = median_obs),
    colour = "#2171b5", size = 0.9, alpha = 0.75
  ) +
  
  ## -- Vertical line: end of observed data
  geom_vline(
    xintercept = M6_DAY,
    linetype = "dotted", colour = "#2171b5", linewidth = 0.7
  ) +
  annotate("text",
           x = M6_DAY + 8, y = y_max * 0.68,
           label = "Last observed\nconcentration (M6)",
           hjust = 0, size = 2.8, colour = "#2171b5") +
  
  ## -- Simulated trough curve (M0–M27)
  geom_line(
    data = pk_pred_daily,
    aes(x = day_round, y = pred),
    colour = "#d95f02", linewidth = 0.75, alpha = 0.90
  ) +
  
  ## -- Clinical timepoint vertical lines
  geom_vline(
    xintercept = c(0, M3_DAY, M27_DAY),
    linetype = "dashed", colour = "grey40", linewidth = 0.5
  ) +
  annotate("text",
           x      = c(2, M3_DAY + 2, M27_DAY - 8),
           y      = y_max,
           vjust  = 1.5,
           hjust  = 0,
           label  = c("M0", "M3", "M27"),
           size   = 3.5, colour = "grey25", fontface = "italic") +
  
  ## -- Cumulative AUROC annotations
  annotate("text",
           x     = M3_DAY / 2,
           y     = y_max * 0.93,
           label = sprintf("Cumulative AUROC\n0→M3\n%.0f ng/mL·day", auroc_M3),
           size  = 2.8, colour = "#2171b5",
           fontface = "bold", hjust = 0.5) +
  annotate("text",
           x     = (M6_DAY + M27_DAY) / 2,
           y     = y_max * 0.93,
           label = sprintf("Cumulative AUROC\n0→M27\n%.0f ng/mL·day", auroc_M27),
           size  = 2.8, colour = "#d95f02",
           fontface = "bold", hjust = 0.5) +
  
  ## -- Scales
  coord_cartesian(ylim = c(0, y_max)) +
  scale_x_continuous(
    name   = "Time (days)",
    breaks = c(0, 90, 180, 270, 360, 450, 540, 630, 720, 810),
    labels = c("0\n(M0)", "90\n(M3)", "180\n(M6)",
               "270", "360", "450",
               "540", "630", "720", "810\n(M27)")
  ) +
  scale_y_continuous(
    name   = "Ibrutinib concentration (ng/mL)",
    expand = expansion(mult = c(0, 0.02))
  ) +
  
  ## -- Labels
  labs(
    title    = "Ibrutinib PK profile — observed vs simulated",
    subtitle = paste0(
      "Blue points/ribbon: observed daily median \u00b1 IQR (M0\u2013M6 only).  ",
      "Red line: simulated daily trough (M0\u2013M27).\n",
      "Shaded areas: cumulative AUROC integration windows ",
      "(blue = 0\u2192M3; orange = 0\u2192M27)."
    )
  ) +
  
  ## -- Theme
  theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 9, colour = "grey40"),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(size = 8)
  )

## ============================================================
## 5. SAVE
## ============================================================

ggsave(file.path(OUTPUT_DIR, "pk_profile_M6.pdf"),
       p_pk, width = 11, height = 5, device = cairo_pdf)

ggsave(file.path(OUTPUT_DIR, "pk_profile_M6.png"),
       p_pk, width = 11, height = 5, dpi = 300)

message("\nSaved:")
message("  ", file.path(OUTPUT_DIR, "pk_profile.pdf"))
message("  ", file.path(OUTPUT_DIR, "pk_profile.png"))

message("\n===== AUROC VALUES FOR ML PIPELINE =====")
message(sprintf("  AUROC_M3  <- %.1f   # cumulative 0\u2192M3  (ng/mL\u00b7day)", auroc_M3))
message(sprintf("  AUROC_M27 <- %.1f  # cumulative 0\u2192M27 (ng/mL\u00b7day)", auroc_M27))
message("=========================================")