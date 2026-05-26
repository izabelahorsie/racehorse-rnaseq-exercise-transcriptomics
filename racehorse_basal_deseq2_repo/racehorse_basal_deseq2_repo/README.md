# Basal/control-only DESeq2 analysis of racehorse whole-blood RNA-seq data

This repository contains the R scripts used to perform basal/control-only differential expression analysis for whole-blood RNA-seq data from Arabian and Thoroughbred racehorses.

The code preserves the original analysis logic and replaces local absolute paths with relative paths defined in configuration files.

## Analyses included

1. **Breed comparison in control-only samples**
   - T1: Thoroughbred vs Arabian
   - T2: Thoroughbred vs Arabian
   - R: Thoroughbred vs Arabian

2. **Performance comparison in control-only race samples**
   - R: podium vs bottom finishers

## Repository structure

```text
racehorse_basal_deseq2_repo/
├── README.md
├── config/
│   ├── config_basal_breed.csv
│   └── config_basal_performance.csv
├── scripts/
│   ├── 00_install_dependencies.R
│   ├── 01_helpers.R
│   ├── 02_plot_functions.R
│   └── 03_run_basal_control_DE.R
├── data/
│   ├── raw/
│   └── metadata/
├── figures/
│   └── example_outputs/
└── results/
```

## Required input files

Place raw count matrices in:

```text
data/raw/
```

Expected files:

```text
rawcounts-matrix_R.csv
rawcounts-matrix_T1.csv
rawcounts-matrix_T2.csv
```

Place metadata/design files in:

```text
data/metadata/
```

Expected files:

```text
design_R_balanced_Arabian_Thoroughbred.csv
design_T1_balanced_Arabian_Thoroughbred.csv
design_T2_balanced_Arabian_Thoroughbred.csv
design_R_with_HorseID_with_performance.xlsx
```

The exact paths are defined in:

```text
config/config_basal_breed.csv
config/config_basal_performance.csv
```

If your filenames differ, update the corresponding config file rather than editing the analysis script.

## How to run

From the repository root, run:

```r
source("scripts/03_run_basal_control_DE.R")
```

The outputs will be written to:

```text
results/Basal_CONTROL_Breed_DESeq2/
results/Basal_CONTROL_Performance_DESeq2/
```

Each analysis folder contains:

- `DE_FULL.csv`
- `DE_SIG_log2FC1_FDR0.05.csv`
- `TOP30_byFDR.csv`
- `PCA.png`
- `Volcano.png`
- `Heatmap_Top30.png`

## Software requirements

The scripts require:

- R
- DESeq2
- dplyr
- stringr
- ggplot2
- pheatmap
- readxl

An optional setup script is provided:

```r
source("scripts/00_install_dependencies.R")
```

## Notes on reproducibility

The workflow performs:

- matching of count matrix columns to metadata samples,
- normalization of condition, breed, and performance labels,
- filtering to control-only samples,
- filtering of genes with total counts below 10,
- DESeq2 differential expression analysis,
- PCA visualization from VST values,
- capped volcano plots,
- heatmaps of top 30 genes ranked by adjusted p-value.

No local computer-specific paths are used in this repository version.
