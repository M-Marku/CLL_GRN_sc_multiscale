# CLL ibrutinib resistance — R pipeline (reorganized)

## Structure

```
CLL_pipeline_R/
├── config.R                          # all paths, colours, shared constants — edit here, not per-script
├── functions.R                       # functions used in >1 script (AUCell cleanup, manual RSS, gprofiler_plot)
└── R/
    ├── 01_extract_CLL_compartment.R  # was: seurat_manipulation_DataExtract.R
    ├── 02_visualise_regulons.R       # was: AUCell_visualisation.R
    ├── 03_RSS_analysis.R             # was: RSS_calculation.R
    ├── 04_fate_regulon_correlation.R # was: corr_aucell_fateProb.R
    └── 05_clonal_phylogenetics.R     # was: cll_phylo_pipeline.R
```

Run in order: each numbered script's output feeds the next. `01` produces the
CLL-only object and the coding-gene matrix pySCENIC needs (pySCENIC itself is
a separate Python step, not included here — send that over when you're ready
and I'll fold it into this same structure). `02`–`04` all consume the pySCENIC
AUCell output. `05` is self-contained (only needs the CLL-only Seurat object +
the full PBMC object) and produces manuscript Figure 4 a–e.

Every script starts with `source("config.R"); source("functions.R")` — run
from a working directory where those two files are visible (i.e., treat
`CLL_pipeline_R/` as your project root, or adjust the `source()` paths).

## What changed, beyond moving files

**Fixed — please check before rerunning downstream analyses:**
- `01_extract_CLL_compartment.R`: `cell_types == c("CLL","Bcells")` →
  `cell_types %in% c("CLL","Bcells")`. The `==` version silently mis-selects
  cells (element-wise recycling against a length-2 vector) rather than
  selecting everything matching either label. This is the very first
  filtering step in the whole pipeline — if this bug was live in results
  you've already generated, everything downstream was built on a subtly
  wrong CLL compartment. Worth rerunning from here and comparing.

**Consolidated — was fragmented/non-running in the original:**
- `03_RSS_analysis.R`: the original `RSS_calculation.R` had several
  redefinitions of `rss_df` and referenced undefined objects
  (`rss_2groups_df`, `label_TFs`, `rss_M27_df`, `rss_vec`, `groups`) —
  looked like leftover fragments from earlier iterations, and wouldn't run
  top-to-bottom. Consolidated into one flow matching what looks like your
  final intent. Flag me if I dropped something you still need from one of
  the earlier versions.

**Added — functionality you described as planned but not yet in the scripts:**
- `03_RSS_analysis.R`, Section B: ranked-based GSEA (fgsea, preranked on RSS)
  as an alternative to the top-10-TF + gprofiler ORA, per your request to
  avoid the arbitrary cutoff. **Both versions are kept and run
  independently** — Section A (unchanged, same figure filename/size as
  before) is what's currently in the manuscript; Section B is new. Nothing
  is silently replaced — decide which one you want in the final figure.
- `04_fate_regulon_correlation.R`, Section B: pathway enrichment of the
  fate-associated regulons (M3/M27), using `gprofiler_plot()` as you
  described — this call didn't exist in the original script, only the
  function definition did (in `AUCell_visualisation.R`).
- `functions.R`: `calc_rss_manual()` — replaces `SCENIC::calcRSS()`, which is
  no longer installable in R. Same Jensen-Shannon-divergence-based formula.
  **Not cross-validated against the original package's output** — if you
  have any old RSS results computed before SCENIC broke, it's worth a quick
  sanity check that the top-regulon ranking matches.

**Standardised:**
- One path was `../Plots_Results/RSS_plot_M27.svg`, inconsistent with every
  other script's `Plots/CLL_only/...` — moved under the same `PLOTS_DIR`.
  Figure dimensions themselves are unchanged (5 × 7).
- Duplicated logic (AUCell column-name cleanup, timepoint extraction from
  barcode) that was copy-pasted near-identically across 3 scripts is now
  `clean_aucell_colnames()` / `extract_timepoint()` in `functions.R`.

## All existing figure dimensions are unchanged

Every `ggsave()`/`svg()` call that existed in your original scripts keeps its
exact `width`/`height`/`dpi` — I did not touch any of these, per your note
that the manuscript figures are already built around them. The only new plots
(ranked GSEA, fate-regulon enrichment) have no prior size to preserve since
they didn't exist before; I picked reasonable defaults and flagged them
clearly in-code as new/adjustable.

## Still to fold in

You mentioned you'll send the rest later — once you do, this is where they'll
plug in:
- pySCENIC / GRNBoost2 / RcisTarget / AUCell pipeline (Python) — produces the
  `AUcell_CLL_only.csv` that scripts 02–04 consume
- CellRank trajectory inference (Python) — produces
  `cell_matrix_fate_probabilities_light.csv` for script 04
- The MaBoSS/PK-coupled notebook (already reviewed separately)

Given the mix of R and Python, the eventual full structure will probably want
a top-level split — something like `R/` and `python/`, each with their own
numbered stages and their own environment file — happy to fold everything
into one coherent tree (with a single top-level README describing the whole
pipeline end-to-end) once I've seen the rest.
