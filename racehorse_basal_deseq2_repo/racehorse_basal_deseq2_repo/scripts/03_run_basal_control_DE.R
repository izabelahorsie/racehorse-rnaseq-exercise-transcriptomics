# Basal/control-only differential expression analysis for racehorse whole-blood RNA-seq data.
# This script preserves the original analysis logic and replaces local absolute paths with config files.

source("scripts/01_helpers.R")
source("scripts/02_plot_functions.R")

# ----------------- CORE DE RUNNER (CONTROL-ONLY) -----------------

run_basal_control_DE <- function(prefix,
                                 counts_path,
                                 design_path,
                                 design_type = c("csv", "xlsx"),
                                 comparison = c("breed", "performance"),
                                 out_base_dir,
                                 lfc_sig = 1,
                                 padj_sig = 0.05){
  
  design_type <- match.arg(design_type)
  comparison <- match.arg(comparison)
  
  message("\n--- Basal CONTROL-only DE: ", prefix, " | ", comparison, " ---")
  
  counts <- read.csv(counts_path, row.names=1, check.names=FALSE)
  
  if (design_type == "csv") {
    meta <- read.csv(design_path, check.names=FALSE)
  } else {
    meta <- readxl::read_xlsx(design_path) %>% as.data.frame(check.names=FALSE)
  }
  
  # Detect common columns
  sample_col <- find_col(meta, c("Sample","SampleID","sample_id","SampleName","Sample_Code","SampleCode","RNA_ID"))
  
  # Condition (for filtering to CONTROL only)
  cond_col <- find_col(meta, c("Condition","condition","Group","CaseControl","case_control"))
  
  # Optional Horse_ID
  horse_col <- NULL
  try(horse_col <- find_col(meta, c("Horse_ID","HorseID","horse_id","horseID","ID_Horse","Animal_ID")), silent=TRUE)
  
  # Breed/performance columns depending on requested comparison
  if (comparison == "breed") {
    breed_col <- find_col(meta, c("Breed","breed"))
    keep_cols <- c(sample_col, cond_col, breed_col)
    if (!is.null(horse_col)) keep_cols <- c(keep_cols, horse_col)
    meta <- meta[, keep_cols, drop=FALSE]
    colnames(meta) <- c("Sample","Condition","Breed", if (!is.null(horse_col)) "Horse_ID" else NULL)
  } else {
    perf_col <- find_col(meta, c("Performance","performance","RankingGroup","ResultGroup","PodiumBottom"))
    keep_cols <- c(sample_col, cond_col, perf_col)
    if (!is.null(horse_col)) keep_cols <- c(keep_cols, horse_col)
    meta <- meta[, keep_cols, drop=FALSE]
    colnames(meta) <- c("Sample","Condition","Performance", if (!is.null(horse_col)) "Horse_ID" else NULL)
  }
  
  meta$Sample <- as.character(meta$Sample)
  
  # Intersect + reorder
  meta <- meta[meta$Sample %in% colnames(counts), , drop=FALSE]
  counts <- counts[, meta$Sample, drop=FALSE]
  meta <- meta[match(colnames(counts), meta$Sample), , drop=FALSE]
  rownames(meta) <- meta$Sample
  
  # Normalize Condition and filter to CONTROL ONLY
  meta$Condition <- factor(normalize_cond(meta$Condition), levels=c("control","case"))
  meta <- meta[meta$Condition == "control", , drop=FALSE]
  counts <- counts[, rownames(meta), drop=FALSE]
  
  if (nrow(meta) < 4) stop("Too few CONTROL samples after filtering.")
  if (comparison == "breed") {
    meta$Breed <- factor(normalize_breed(meta$Breed), levels=c("arabian","thoroughbred"))
    if (length(unique(meta$Breed)) < 2) stop("Need both breeds in CONTROL samples.")
    group_col <- "Breed"
    left_level <- "arabian"
    right_level <- "thoroughbred"
    design_formula <- ~ Breed
    contrast <- c("Breed", "thoroughbred", "arabian")  # positive = TB higher than Arabian
  } else {
    meta$Performance <- factor(normalize_perf(meta$Performance), levels=c("bottom","podium"))
    if (length(unique(meta$Performance)) < 2) stop("Need both performance groups in CONTROL samples.")
    group_col <- "Performance"
    left_level <- "podium"
    right_level <- "bottom"
    design_formula <- ~ Performance
    contrast <- c("Performance", "podium", "bottom")   # positive = podium higher than bottom
  }
  
  # Basic count filter
  keep <- rowSums(counts) >= 10
  counts <- counts[keep, , drop=FALSE]
  
  # DESeq2
  dds <- DESeqDataSetFromMatrix(countData=round(as.matrix(counts)),
                                colData=meta,
                                design=design_formula)
  dds <- DESeq(dds)
  
  res <- results(dds, contrast=contrast)
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df <- res_df %>% arrange(padj)
  
  sig_df <- res_df %>%
    filter(!is.na(padj), padj < padj_sig, abs(log2FoldChange) > lfc_sig)
  
  top30 <- res_df %>%
    filter(!is.na(padj)) %>%
    slice_head(n=30) %>%
    pull(gene)
  
  # Output folder
  out_dir <- file.path(out_base_dir, paste0(prefix, "_CONTROL_", comparison))
  make_dir(out_dir)
  
  # Save tables
  write.csv(res_df, file.path(out_dir, "DE_FULL.csv"), row.names=FALSE)
  write.csv(sig_df, file.path(out_dir, paste0("DE_SIG_log2FC", lfc_sig, "_FDR", padj_sig, ".csv")), row.names=FALSE)
  write.csv(data.frame(gene=top30), file.path(out_dir, "TOP30_byFDR.csv"), row.names=FALSE)
  
  # VST for plots
  vsd <- vst(dds, blind=FALSE)
  vst_mat <- assay(vsd)
  
  # PCA
  p_pca <- pca_plot(vst_mat, meta, group_col, paste0(prefix, " | CONTROL | PCA: ", group_col))
  save_gg(p_pca, file.path(out_dir, "PCA.png"), 7, 5)
  
  # Volcano
  p_vol <- volcano_plot_capped(res_df, paste0(prefix, " | CONTROL | Volcano: ", comparison),
                               x_cap=5, y_cap=20)
  save_gg(p_vol, file.path(out_dir, "Volcano.png"), 7, 5)
  
  # Heatmap top30 (gene clustering only; samples fixed order)
  if (length(top30) >= 2) {
    heatmap_top30(vst_mat, meta, top30,
                  out_png = file.path(out_dir, "Heatmap_Top30.png"),
                  title = paste0(prefix, " | CONTROL | Top30 | ", left_level, " (left) vs ", right_level, " (right)"),
                  group_col = group_col,
                  left_level = left_level,
                  right_level = right_level)
  }
  
  message("✅ Saved: ", out_dir)
  invisible(list(full=res_df, sig=sig_df, top30=top30))
}

# ----------------- RUN ALL REQUESTED ANALYSES -----------------

run_config <- function(config_path, results_base_dir = "results"){
  cfg <- read.csv(config_path, stringsAsFactors=FALSE, check.names=FALSE)
  required_cols <- c("dataset", "counts", "design", "design_type", "comparison", "out_subdir")
  missing_cols <- setdiff(required_cols, colnames(cfg))
  if (length(missing_cols) > 0) {
    stop("Missing required config columns: ", paste(missing_cols, collapse=", "))
  }
  
  for (i in seq_len(nrow(cfg))) {
    out_dir <- file.path(results_base_dir, cfg$out_subdir[i])
    make_dir(out_dir)
    
    run_basal_control_DE(
      prefix = cfg$dataset[i],
      counts_path = cfg$counts[i],
      design_path = cfg$design[i],
      design_type = cfg$design_type[i],
      comparison = cfg$comparison[i],
      out_base_dir = out_dir
    )
  }
}

run_config("config/config_basal_breed.csv")
run_config("config/config_basal_performance.csv")

cat("\nAll basal CONTROL-only analyses completed.\n")
