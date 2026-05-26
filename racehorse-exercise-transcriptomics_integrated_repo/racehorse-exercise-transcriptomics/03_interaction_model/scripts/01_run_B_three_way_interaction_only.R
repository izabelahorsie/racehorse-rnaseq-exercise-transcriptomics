suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readxl)
  library(openxlsx)
})

# ----------------- HELPERS -----------------

find_col <- function(df, candidates){
  cn <- colnames(df)
  hit <- cn[tolower(cn) %in% tolower(candidates)]
  if (length(hit)>0) return(hit[1])
  for (cand in candidates){
    hit <- cn[grepl(tolower(cand), tolower(cn))]
    if (length(hit)>0) return(hit[1])
  }
  stop(paste("Missing column among:", paste(candidates, collapse=", ")))
}

normalize_cond <- function(x){
  z <- tolower(trimws(as.character(x)))
  ifelse(z %in% c("control","ctrl") | grepl("con", z), "control", "case")
}

normalize_breed <- function(x){
  z <- tolower(trimws(as.character(x)))
  ifelse(grepl("arab|ara", z), "arabian",
         ifelse(grepl("thor|thorough", z), "thoroughbred", z))
}

normalize_perf <- function(x){
  z <- tolower(trimws(as.character(x)))
  ifelse(z %in% c("podium","top","winner","high","best") | grepl("pod", z), "podium",
         ifelse(z %in% c("bottom","low","worst","poor") | grepl("bot", z), "bottom", z))
}

make_dir <- function(path){
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

res_table <- function(dds, name){
  res <- results(dds, name = name)
  df <- as.data.frame(res)
  df$gene <- rownames(df)
  df <- df %>% arrange(padj)
  df
}

# ----------------- CORE FUNCTION -----------------

run_B_three_way_interaction_only <- function(counts_path,
                                             design_xlsx,
                                             out_dir,
                                             dataset = "R"){
  make_dir(out_dir)
  
  counts <- read.csv(counts_path, row.names = 1, check.names = FALSE)
  meta <- readxl::read_xlsx(design_xlsx) %>% as.data.frame(check.names = FALSE)
  
  sample_col <- find_col(meta, c("Sample","SampleID","sample_id","SampleName","Sample_Code","SampleCode","RNA_ID"))
  cond_col   <- find_col(meta, c("Condition","condition","Group","CaseControl","case_control"))
  breed_col  <- find_col(meta, c("Breed","breed"))
  perf_col   <- find_col(meta, c("Performance","performance","RankingGroup","ResultGroup","PodiumBottom"))
  horse_col  <- find_col(meta, c("Horse_ID","HorseID","horse_id","horseID","ID_Horse","Animal_ID"))
  
  meta <- meta[, c(sample_col, cond_col, breed_col, perf_col, horse_col)]
  colnames(meta) <- c("Sample","Condition","Breed","Performance","Horse_ID")
  meta$Sample <- as.character(meta$Sample)
  
  # Intersect and align metadata and count matrix.
  meta <- meta[meta$Sample %in% colnames(counts), , drop = FALSE]
  counts <- counts[, meta$Sample, drop = FALSE]
  meta <- meta[match(colnames(counts), meta$Sample), , drop = FALSE]
  rownames(meta) <- meta$Sample
  
  # Normalize factor labels and preserve the original reference levels.
  meta$Condition   <- factor(normalize_cond(meta$Condition), levels = c("control","case"))
  meta$Breed       <- factor(normalize_breed(meta$Breed), levels = c("arabian","thoroughbred"))
  meta$Performance <- factor(normalize_perf(meta$Performance), levels = c("bottom","podium"))
  meta$Horse_ID    <- factor(meta$Horse_ID)
  
  if (length(unique(meta$Breed)) < 2) stop("Need both breeds present.")
  if (length(unique(meta$Performance)) < 2) stop("Need both performance groups present.")
  if (length(unique(meta$Condition)) < 2) stop("Need both conditions present.")
  
  # Keep the original filtering rule.
  keep <- rowSums(counts) >= 10
  counts <- counts[keep, , drop = FALSE]
  
  design_formula <- ~ Breed * Performance * Condition
  message("Using design: ", deparse(design_formula))
  
  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(counts)),
    colData = meta,
    design = design_formula
  )
  
  dds <- DESeq(dds, test = "Wald")
  
  rn <- resultsNames(dds)
  message("Model coefficients:\n", paste(rn, collapse = "\n"))
  
  int_BPC <- rn[
    grepl("Breed.*Performance.*Condition|Breed.*Condition.*Performance|Performance.*Breed.*Condition|Performance.*Condition.*Breed|Condition.*Breed.*Performance|Condition.*Performance.*Breed", rn)
  ][1]
  
  if (is.na(int_BPC)) {
    stop("Three-way interaction term not found in resultsNames(dds). Check factor levels and model coefficients.")
  }
  
  df_int <- res_table(dds, name = int_BPC)
  
  wb <- createWorkbook()
  addWorksheet(wb, "INTERACTION_Breed_x_Perf_x_Cond")
  writeData(wb, "INTERACTION_Breed_x_Perf_x_Cond", df_int)
  
  out_xlsx <- file.path(
    out_dir,
    "DESeq2_B_Breed_x_Performance_interaction_results.xlsx"
  )
  
  saveWorkbook(wb, out_xlsx, overwrite = TRUE)
  
  # Lightweight reproducibility records.
  writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
  write.csv(
    data.frame(
      dataset = dataset,
      design = deparse(design_formula),
      coefficient_exported = int_BPC,
      n_samples = ncol(counts),
      n_genes_after_filtering = nrow(counts),
      output_file = out_xlsx
    ),
    file.path(out_dir, "analysis_summary.csv"),
    row.names = FALSE
  )
  
  message("Saved only three-way interaction sheet: ", out_xlsx)
  
  invisible(list(
    dds = dds,
    interaction_name = int_BPC,
    results = df_int,
    out_xlsx = out_xlsx
  ))
}

# ----------------- RUN FROM CONFIG -----------------

config_path <- "config/config_interaction_B_only.csv"
config <- read.csv(config_path, stringsAsFactors = FALSE)

for (i in seq_len(nrow(config))) {
  run_B_three_way_interaction_only(
    counts_path = config$counts[i],
    design_xlsx = config$design[i],
    out_dir = config$out_dir[i],
    dataset = config$dataset[i]
  )
}

cat("\nAll three-way interaction analyses completed.\n")
