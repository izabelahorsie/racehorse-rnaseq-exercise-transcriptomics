# Racehorse exercise transcriptomics

Reproducible RNA-seq analysis workflows and processed count data structure for investigating exercise-, breed-, and performance-associated whole-blood transcriptomic responses in Arabian and Thoroughbred racehorses.

## Repository structure

```text
racehorse-exercise-transcriptomics/
├── shared_data/
│   ├── raw/
│   └── metadata/
├── 01_breed_stratified_DE/
├── 02_basal_control_DE/
├── 03_interaction_model/
└── 04_performance_stratified_DE/
```

## Required input files

Place count matrices in:

```text
shared_data/raw/
```

Expected files:

```text
rawcounts-matrix_R.csv
rawcounts-matrix_T1.csv
rawcounts-matrix_T2.csv
```

Place design/metadata files in:

```text
shared_data/metadata/
```

Expected files:

```text
design_R_balanced_Arabian_Thoroughbred.csv
design_T1_balanced_Arabian_Thoroughbred.csv
design_T2_balanced_Arabian_Thoroughbred.csv
design_R_with_HorseID_with_performance.xlsx
```

## Analyses

### 01_breed_stratified_DE
Breed-stratified DESeq2 analysis of exercise response within Arabian and Thoroughbred horses across R, T1 and T2 datasets.

Run from inside the folder:

```r
setwd("01_breed_stratified_DE")
source("scripts/02_run_breed_deseq2.R")
```

### 02_basal_control_DE
Control-only basal transcriptomic comparisons for breed and performance analyses.

Run from inside the folder:

```r
setwd("02_basal_control_DE")
source("scripts/03_run_basal_control_DE.R")
```

### 03_interaction_model
Three-way Breed × Performance × Condition interaction model. The script exports only the `INTERACTION_Breed_x_Perf_x_Cond` sheet in `DESeq2_B_Breed_x_Performance_interaction_results.xlsx`.

Run from inside the folder:

```r
setwd("03_interaction_model")
source("scripts/01_run_B_three_way_interaction_only.R")
```

### 04_performance_stratified_DE
Performance-stratified DESeq2 analysis of post-race exercise response within podium and bottom finishers.

Run from inside the folder:

```r
setwd("04_performance_stratified_DE")
source("scripts/03_run_performance_stratified_condition_DE.R")
```

## R package requirements

The scripts use:

- DESeq2
- dplyr
- stringr
- ggplot2
- pheatmap
- readxl
- openxlsx

A reproducible `renv.lock` can be generated after installing the required packages:

```r
install.packages("renv")
renv::init()
renv::snapshot()
```

## Notes

The repository is organized with one shared input-data location to avoid duplicating count matrices and metadata across analysis folders. The scripts preserve the original analytical logic and replace local Windows paths with relative repository paths defined in each analysis-specific `config/` file.
