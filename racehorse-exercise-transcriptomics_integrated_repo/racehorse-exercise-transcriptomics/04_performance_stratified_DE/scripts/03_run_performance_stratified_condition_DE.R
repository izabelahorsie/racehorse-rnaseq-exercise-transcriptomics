suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readxl)
})

source("scripts/01_helpers.R")
source("scripts/02_plot_functions.R")

# choose design: use pairing if Horse_ID has both control and case
choose_design_condition <- function(meta){
  tab <- table(meta$Horse_ID, meta$Condition)
  has_pairs <- any(rowSums(tab > 0) == 2)
  if (has_pairs) as.formula("~ Horse_ID + Condition") else as.formula("~ Condition")
}

# ----------------- CORE: RUN DESeq2 WITHIN ONE PERFORMANCE GROUP -----------------

run_deseq2_within_performance <- function(counts, meta, perf_level, out_dir){
  make_dir(out_dir)
  
  sub_meta <- meta %>% filter(Performance == perf_level)
  if (nrow(sub_meta) < 4) {
    message("Skipping ", perf_level, ": too few samples.")
    return(invisible(NULL))
  }
  if (length(unique(sub_meta$Condition)) < 2) {
    message("Skipping ", perf_level, ": only one condition present.")
    return(invisible(NULL))
  }
  
  sub_counts <- counts[, rownames(sub_meta), drop=FALSE]
  
  # basic count filter
  keep <- rowSums(sub_counts) >= 10
  sub_counts <- sub_counts[keep, , drop=FALSE]
  
  design_formula <- choose_design_condition(sub_meta)
  message("Performance=", perf_level, " | Using design: ", deparse(design_formula))
  
  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(sub_counts)),
    colData   = sub_meta,
    design    = design_formula
  )
  dds$Horse_ID  <- droplevels(dds$Horse_ID)
  dds$Condition <- droplevels(dds$Condition)
  
  dds <- DESeq(dds)
  
  res <- results(dds, contrast=c("Condition","case","control"))
  
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df <- res_df %>% arrange(padj)
  
  sig_df <- res_df %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)
  
  top30 <- res_df %>%
    filter(!is.na(padj)) %>%
    slice_head(n=30) %>%
    pull(gene)
  
  # Save tables
  write.csv(res_df, file.path(out_dir, paste0("DESeq2_", perf_level, "_FULL_control_vs_case.csv")), row.names=FALSE)
  write.csv(sig_df, file.path(out_dir, paste0("DESeq2_", perf_level, "_SIG_log2FC1_FDR0.05_control_vs_case.csv")), row.names=FALSE)
  write.csv(data.frame(gene=top30),
            file.path(out_dir, paste0("DESeq2_", perf_level, "_TOP30_byFDR_control_vs_case.csv")),
            row.names=FALSE)
  
  # VST for plots
  vsd <- vst(dds, blind=FALSE)
  vst_mat <- assay(vsd)
  
  # PCA
  p_pca <- pca_plot(vst_mat, sub_meta, paste0("R | ", perf_level, " | PCA (control vs case)"))
  save_gg(p_pca, file.path(out_dir, paste0("PCA_", perf_level, ".png")))
  
  # Volcano
  p_vol <- volcano_plot_capped(res_df, paste0("R | ", perf_level, " | Volcano (control vs case)"),
                               x_cap=5, y_cap=60)
  save_gg(p_vol, file.path(out_dir, paste0("Volcano_", perf_level, ".png")))
  
  # Heatmap top30
  if (length(top30) >= 2) {
    heatmap_top30(vst_mat, sub_meta, top30,
                  out_png = file.path(out_dir, paste0("Heatmap_Top30_", perf_level, "_caseLEFT_controlRIGHT.png")),
                  title   = paste0("R | ", perf_level, " | Top30 (case left, control right; gene clustering only)"))
  }
  
  invisible(list(full=res_df, sig=sig_df, top30=top30))
}

# ----------------- MAIN WORKFLOW -----------------

run_performance_stratified_condition_DE <- function(counts_path, design_xlsx, out_base_dir){
  
  make_dir(out_base_dir)
  
  # read counts
  counts <- read.csv(counts_path, row.names=1, check.names=FALSE)
  
  # read design from xlsx
  meta <- readxl::read_xlsx(design_xlsx) %>% as.data.frame(check.names=FALSE)
  
  # detect columns
  sample_col <- find_col(meta, c("Sample","SampleID","sample_id","SampleName","Sample_Code","SampleCode","RNA_ID"))
  cond_col   <- find_col(meta, c("Condition","condition","Group","CaseControl","case_control"))
  perf_col   <- find_col(meta, c("Performance","performance","RankingGroup","ResultGroup","PodiumBottom"))
  horse_col  <- find_col(meta, c("Horse_ID","HorseID","horse_id","horseID","ID_Horse","Animal_ID"))
  
  meta <- meta[, c(sample_col, cond_col, perf_col, horse_col)]
  colnames(meta) <- c("Sample","Condition","Performance","Horse_ID")
  meta$Sample <- as.character(meta$Sample)
  
  # intersect + align
  meta <- meta[meta$Sample %in% colnames(counts), , drop=FALSE]
  counts <- counts[, meta$Sample, drop=FALSE]
  meta <- meta[match(colnames(counts), meta$Sample), , drop=FALSE]
  rownames(meta) <- meta$Sample
  
  # normalize labels
  meta$Condition   <- factor(normalize_cond(meta$Condition), levels=c("control","case"))
  meta$Performance <- factor(normalize_perf(meta$Performance), levels=c("bottom","podium"))
  meta$Horse_ID    <- factor(meta$Horse_ID)
  
  # run for each performance group separately
  for (perf in c("podium","bottom")){
    out_dir <- file.path(out_base_dir, paste0("R_", perf, "_control_vs_case"))
    message("\n--- Running subgroup: ", perf, " ---")
    run_deseq2_within_performance(counts, meta, perf, out_dir)
  }
  
  message("\n✅ Done. Output base dir: ", out_base_dir)
}

# ----------------- CONFIG-BASED RUN -----------------

config <- read.csv("config/config_performance_stratified.csv", stringsAsFactors = FALSE)

for (i in seq_len(nrow(config))) {
  run_performance_stratified_condition_DE(
    counts_path = config$counts[i],
    design_xlsx = config$design[i],
    out_base_dir = config$output_dir[i]
  )
}
