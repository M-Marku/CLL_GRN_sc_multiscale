# A multiscale pharmacokinetic Boolean framework for ibrutinib resistance in a CLL patient case study

**Malvina Marku, Cécile Pineau, Marcelo Hurtado, Hugo Chenel, Vera Pancaldi, Jean-Baptiste Woillard, Loïc Ysebaert, Étienne Chatelut, Melanie White-Koning**

*Centre de Recherches en Cancérologie de Toulouse (CRCT)*

## Overview

Code and model files for a case study that couples a patient-specific ibrutinib pharmacokinetic (PK) model with single-cell gene regulatory network (GRN) inference and a stochastic Boolean model of CLL signalling, using the longitudinal CITE-seq data of one patient (Cadot et al., 2020).

Analyses: GRN inference (pySCENIC) · trajectory inference (CellRank) · PK-coupled Boolean model (MaBoSS) · PK–expression association.

## Data and results

Processed data and results are on Zenodo: [Zenodo DOI] (`data.zip` and `results.zip`). Extract them into `data/` and `results/`.

- **Original CITE-seq data:** Cadot et al. (2020), *Biomarker Research* 8:72, https://doi.org/10.1186/s40364-020-00253-w
- **Patient PK data:** not publicly shareable; contact malvina.marku@inserm.fr.
- **pySCENIC reference files:** downloaded from https://resources.aertslab.org/cistarget/ (`hg19-500bp-upstream-7species.mc9nr.feather`, `hg19-tss-centered-5kb-7species.mc9nr.feather`, `motifs-v9-nr.hgnc-m0.001-o0.0.tbl`, `allTFs_hg38.txt`).

## Repository structure

```
├── model/       Boolean model (MaBoSS .bnd, .cfg, .graphml)
├── scripts/
│   ├── notebooks/   pySCENIC, CellRank, PK-Boolean simulation
│   └── R/           data extraction, RSS and centrality, AUCell binarisation,
│                    Slingshot, PK profile, PK–expression association
├── data/        see "Data and results"
└── results/     outputs of the scripts (on Zenodo)
```

## Running the analyses

Run in this order:

1. `scripts/R/01_seurat_data_extraction.R`
2. `scripts/notebooks/01_GRN_pySCENIC.ipynb`
3. `scripts/R/02_RSS_calculation.R`, `03_centrality_measures.R`, `05_infer_interaction_signs.R`, `06_AUCell_binarisation.R`
4. `scripts/notebooks/02_trajectory_CellRank.ipynb`, then `scripts/R/07_slingshot_trajectory.R`, `08_corr_aucell_fateprob.R`
5. `scripts/R/10_PK_profile.R`
6. `scripts/notebooks/03_Boolean_PK_coupled.ipynb`
7. `scripts/R/12_grn_pk_correlation.R`

## Requirements

- **Python:** 3.12, run in the [CoLoMoTo Docker](https://github.com/colomoto/colomoto-docker) image (pySCENIC 0.12.1, CellRank 2.0.7, Scanpy 1.10.3, pyMaBoSS). Environment file: `environment/environment_python.yml`.
- **R:** 4.3.x (Seurat, AUCell, igraph, tidyverse, pheatmap, patchwork). Install with `source("environment/requirements_R.R")`.

pySCENIC was run on the Genotoul cluster (https://bioinfo.genotoul.fr); the other analyses run on a standard workstation.

## Citation

Marku M, Pineau C, Hurtado M, Chenel H, Pancaldi V, Woillard J-B, Ysebaert L, Chatelut E, White-Koning M. ([year]). A multiscale pharmacokinetic Boolean framework for ibrutinib resistance in a CLL patient case study. *[Journal]*. doi: [to be added]

Please also cite Cadot et al. (2020).

## Funding, contact and licence

Funded by the PEPR Santé Numérique, ANR-23-CE45-0014, EUR CARe ANR-18-EURE-0003 and a CIFRE fellowship (No. 231967A10, with Evotec France SAS).
Contact: Malvina Marku, malvina.marku@inserm.fr. Licence: MIT (see `LICENSE`).
