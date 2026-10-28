# A multiscale pharmacokinetic Boolean framework for ibrutinib resistance in a CLL patient case study

**Malvina Marku, Cécile Pineau, Marcelo Hurtado, Hugo Chenel, Jean-Baptiste Woillard, Loïc Ysebaert, Étienne Chatelut, Vera Pancaldi, Melanie White-Koning**


## Overview

Code and model files for a case study that couples a patient-specific ibrutinib pharmacokinetic (PK) model with single-cell gene regulatory network (GRN) inference and a stochastic Boolean model of CLL signalling, using the longitudinal CITE-seq data of one patient (Cadot et al., 2020).

Analyses: GRN inference (pySCENIC) · trajectory inference (CellRank) · PK-coupled Boolean model (MaBoSS, Colomoto docker) · PK-expression association.

## Data and results

Processed data and results are on Zenodo: [https://doi.org/10.5281/zenodo.23320593] (`Data.zip` and `Results.zip`). Extract them into `Data/` and `Results/`.

- **Original CITE-seq data:** Cadot et al. (2020), *Biomarker Research* 8:72, https://doi.org/10.1186/s40364-020-00253-w
- **Patient PK data:** not publicly shareable; contact malvina.marku@inserm.fr.
- **pySCENIC reference files:** downloaded from https://resources.aertslab.org/cistarget/ (`hg19-500bp-upstream-7species.mc9nr.feather`, `hg19-tss-centered-5kb-7species.mc9nr.feather`, `motifs-v9-nr.hgnc-m0.001-o0.0.tbl`, `allTFs_hg38.txt`).

## Repository structure

```
├── Model/       Boolean model (MaBoSS .bnd, .cfg)
├── Scripts/
│   ├── notebooks/   pySCENIC, CellRank, PK-Boolean simulation
│   └── R/           data extraction, RSS and centrality, AUCell binarisation,
│                    Slingshot, PK profile, PK–expression association
├── Data/        see "Data and results"
└── Results/     outputs of the scripts
```

## Running the analyses

Run in this order:

1. `scripts/R/01_seurat_data_extraction.R`
2. `scripts/notebooks/pySCENIC - Full pipeline-CLL.ipynb`
3. `scripts/R/03_RSS_calculation.R`
4. `scripts/notebooks/cellRank_trajectory.ipynb`, then `04_fate_regulon_correlation.R`
5. `scripts/notebooks/MaBoSS_PK_coupled_clean.ipynb`
6. `scripts/R/05_grn_pk_correlation_pipeline.R`

## Requirements

- **Python:** 3.12, run in the [CoLoMoTo Docker](https://github.com/colomoto/colomoto-docker) image (pySCENIC 0.12.1, CellRank 2.0.7, Scanpy 1.10.3, pyMaBoSS).
- **R:** 4.3.x (Seurat, AUCell, igraph, tidyverse, pheatmap, patchwork).

pySCENIC was run on the Genotoul cluster (https://bioinfo.genotoul.fr); the other analyses run on a standard workstation.

## Citation

Marku M, Pineau C, Hurtado M, Chenel H, Woillard J-B, Ysebaert L, Chatelut E, Pancaldi V, White-Koning M. (2026). A multiscale pharmacokinetic Boolean framework for ibrutinib resistance in a CLL patient case study. [biorxiv preprint](https://doi.org/10.64898/2026.10.26.759125)

