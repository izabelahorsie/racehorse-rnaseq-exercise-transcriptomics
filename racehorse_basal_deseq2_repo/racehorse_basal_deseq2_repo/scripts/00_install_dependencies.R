# Optional setup script. Run once if required.
# BiocManager is used for DESeq2; CRAN is used for the remaining packages.

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_packages <- c("dplyr", "stringr", "ggplot2", "pheatmap", "readxl")
missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran) > 0) {
  install.packages(missing_cran)
}

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
}
