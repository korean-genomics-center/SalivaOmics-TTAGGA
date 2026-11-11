# SalivaOmics-TTAGGA

Figure and analysis code for:

**Genome-wide DNA methylation concordance between blood and saliva in a single
Korean individual, characterized against a personal telomere-to-telomere genome
assembly: a pilot study**

Choi *et al.* (2026), manuscript under review at *Genomics & Informatics*.

Licensed under the MIT License (see `LICENSE`).

---

## Overview

This repository contains the R scripts that generate the manuscript figures from
the processed methylation tables (per-CpG beta/coverage matrices and DSS DMR
calls). The study profiles blood and saliva from one individual (KOREF1) by
EM-seq and analyses the data against both hg38 and the individual's personal
telomere-to-telomere assembly (KOREF-TTAGGA).

This repository contains **analysis and figure code only**; no sequence data or
intermediate tables are redistributed here.

---

## Data availability

The EM-seq methylation data of blood and saliva generated in this study are
being deposited in the Korea BioData Station (K-BDS), a controlled-access
repository operated by the Korea Bioinformation Center (KOBIC); the accession
number will be provided upon completion. The personal KOREF-TTAGGA assembly is
described in a separate manuscript (currently a preprint under review).

---

## System requirements

**Operating system.** Linux (x86_64); developed on an SGE cluster. The figure
scripts also run on any desktop with R.

**R.** R 4.1 or later. Package versions used for the published figures are listed
in `env/r_packages.txt`: ggplot2 3.5.1, patchwork 1.3.0, data.table 1.16.2,
scales 1.3.0.

```r
install.packages(c("ggplot2", "patchwork", "data.table", "scales"))
```

---

## Instructions for use

### Paths

Figure scripts read their input tables from a single project data root. Set it in
your environment before running, or edit the two lines at the top of each script:

```bash
export SALIVA_ROOT=/path/to/your/dedup_rerun
```

Each script reads processed tables under `${SALIVA_ROOT}` (beta/coverage
matrices, DSS DMR output, enrichment and marker-overlap tables) and writes the
figure PDF/PNG to `${SALIVA_ROOT}/figures_rev`. No data are bundled in this
repository; the scripts are provided so that each published figure can be traced
to the computation that produced it.

### Colours

All figures use a colour-blind-safe palette (Wong 2011, *Nature Methods* 8:441),
defined in `figures/_config.R`. The palette values are those used for the
published figures and should not be changed.

---

## Script-to-figure mapping

Several script file names predate the final figure numbering; the mapping below
is to the figures as numbered in the submitted manuscript.

| Script | Produces |
|---|---|
| `figures/Fig1_QC_rev1.R` | **Figure 1** — sequencing/QC overview, methylation density, coverage thresholds, and blood/saliva CpG sharing |
| `figures/Fig2_DMR_rev1.R` | **Figure 2** — DMR landscape (volcano, length distribution, direction counts, genomic annotation) under hg38 |
| `figures/Fig3_Biology_rev2.R` | **Figure 3** — functional interpretation (GO BP and KEGG on matched background; cell-type-marker overlap, Loyfer atlas) |
| `figures/Fig5_Correction_rev1.R` | **Figure 4** — blood–saliva concordance and held-out correction benchmark (uncorrected, blood-mean null, global linear, global mean-shift) |
| `figures/Fig4_T2T_rev1.R` | **Figure 5** — personal TTAGGA effect (CpG counts, genomic classification under hg38 vs personal annotation, EWAS-Catalog overlap) |
| `figures/FigS_DMR_fate_rev2.R` | **Supplementary Figure 5** — DMR correspondence between references (matched / hg38-only / KOREF-only) |
| `figures/FigB_genomic_classification.R` | Stand-alone genomic-classification panel (Table 2 / Fig. 5B values; `n = 135,098` dedup KOREF-TTAGGA DMRs) |
| `figures/_config.R` | Shared palette and theme imported by the figure scripts |

---

## Contact

- Jong Bhak — jongbhak@genomics.org (corresponding author)
- Hyoungjin Choi — hjchoi349@unist.ac.kr, ORCID 0009-0007-0504-8870 (analysis and figure code)

## Funding

Ulsan Metropolitan City Research Fund (Grant No. 1.200047.01), UNIST; and the
National Research Foundation of Korea (NRF) grants funded by the Korea government
(MSIT) (RS-2026-25603518 and RS-2026-25600210).
