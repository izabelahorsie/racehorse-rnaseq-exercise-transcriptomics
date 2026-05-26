# Breed-stratified DESeq2 analysis runner
# This script preserves the original analytical logic and uses paths from config/config_breed_analysis.csv.

source("scripts/01_helpers.R")

# ----------------- CORE RUNNER -----------------

run_deseq2_dataset_breed <- function(prefix, counts_path, design_path, breed_target, out_base_dir){
  
  message("\n--- DESeq2: Dataset=", prefix, " | Breed=", breed_target, " ---")
  
  counts <- read.csv(counts_path, row.names=1, check.names=FALSE)
  meta   <- read.csv(design_path, check.names=FALSE)
  
  sample_col <- find_col(meta, c("Sample","SampleID","sample_id","SampleName","Sample_Code","SampleCode","RNA_ID"))
  breed_col  <- find_col(meta, c("Breed","breed"))
  cond_col   <- find_col(meta, c("Condition","condition","Group","CaseControl","case_control"))
  horse_col  <- find_col(meta, c("Horse_ID","HorseID","horse_id","horseID","ID_Horse","Animal_ID"))
  
  meta <- meta[, c(sample_col, breed_col, cond_col, horse_col)]
  colnames(meta) <- c("Sample","Breed","Condition","Horse_ID")
  meta$Sample <- as.character(meta$Sample)
  
  # intersect + reorder
  meta <- meta[meta$Sample %in% colnames(counts), , drop=FALSE]
  counts <- counts[, meta$Sample, drop=FALSE]
  meta <- meta[match(colnames(counts), meta$Sample), , drop=FALSE]
  rownames(meta) <- meta$Sample
  
  # normalize labels
  meta$Breed     <- factor(normalize_breed(meta$Breed), levels=c("arabian","thoroughbred"))
  meta$Condition <- factor(normalize_cond(meta$Condition), levels=c("control","case"))
  meta$Horse_ID  <- factor(meta$Horse_ID)
  
  # subset breed
  meta_b <- meta[meta$Breed == breed_target, , drop=FALSE]
  counts_b <- counts[, rownames(meta_b), drop=FALSE]
  
  if (nrow(meta_b) < 4) {
    message("Skipping: too few samples.")
    return(invisible(NULL))
  }
  if (length(unique(meta_b$Condition)) < 2) {
    message("Skipping: only one Condition level.")
    return(invisible(NULL))
  }
  
  # basic count filter
  keep <- rowSums(counts_b) >= 10
  counts_b <- counts_b[keep, , drop=FALSE]
  
  design_formula <- choose_design(meta_b)
  message("Using design: ", deparse(design_formula))
  
  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(counts_b)),
    colData   = meta_b,
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
  
  # output dir
  outdir <- file.path(out_base_dir, paste0(prefix, "_", breed_target))
  make_dir(outdir)
  
  # write tables
  write.csv(res_df, file.path(outdir, paste0("DESeq2_", prefix, "_", breed_target, "_FULL.csv")), row.names=FALSE)
  write.csv(sig_df, file.path(outdir, paste0("DESeq2_", prefix, "_", breed_target, "_SIG_log2FC1_FDR0.05.csv")), row.names=FALSE)
  write.csv(data.frame(gene=top30),
            file.path(outdir, paste0("DESeq2_", prefix, "_", breed_target, "_TOP30_byFDR.csv")),
            row.names=FALSE)
  
  # VST for plots
  vsd <- vst(dds, blind=FALSE)
  vst_mat <- assay(vsd)
  
  # PCA
  p_pca <- pca_plot(vst_mat, meta_b, paste0(prefix, " | ", breed_target, " | PCA (VST)"))
  save_gg(p_pca, file.path(outdir, paste0("PCA_", prefix, "_", breed_target, ".png")), 7, 5)
  
  # Volcano
  p_vol <- volcano_plot(res_df, paste0(prefix, " | ", breed_target, " | Volcano"))
  save_gg(p_vol, file.path(outdir, paste0("Volcano_", prefix, "_", breed_target, ".png")), 7, 5)
  
  # Heatmap top30
  if (length(top30) >= 2) {
    heatmap_top30(
      vst_mat = vst_mat,
      coldata = meta_b,
      top_genes = top30,
      out_png = file.path(outdir, paste0("Heatmap_Top30_", prefix, "_", breed_target, "_caseLEFT_controlRIGHT.png")),
      title = paste0(prefix, " | ", breed_target, " | Top30 (case left, control right; gene clustering only)")
    )
  }
  
  message("✅ Saved to: ", outdir)
  invisible(list(full=res_df, sig=sig_df, top30=top30))
}

# ----------------- INPUT CONFIGURATION -----------------

config_path <- "config/config_breed_analysis.csv"
if (!file.exists(config_path)) {
  stop("Configuration file not found: ", config_path)
}

inputs <- read.csv(config_path, stringsAsFactors = FALSE, check.names = FALSE)
required_cols <- c("dataset", "counts", "design")
missing_cols <- setdiff(required_cols, colnames(inputs))
if (length(missing_cols) > 0) {
  stop("Missing required column(s) in config file: ", paste(missing_cols, collapse = ", "))
}

out_base_dir <- "results/DESeq2_outputs"
make_dir(out_base_dir)

# ----------------- RUN ALL -----------------

for (i in seq_len(nrow(inputs))){
  pref <- inputs$dataset[i]
  for (breed in c("arabian","thoroughbred")){
    run_deseq2_dataset_breed(
      prefix = pref,
      counts_path = inputs$counts[i],
      design_path = inputs$design[i],
      breed_target = breed,
      out_base_dir = out_base_dir
    )
  }
}

cat("\nAll done.\nOutputs in: ", out_base_dir, "\n")
