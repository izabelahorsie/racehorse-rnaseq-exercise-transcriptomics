suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(pheatmap)
})

# ----------------- HELPERS -----------------

normalize_cond <- function(x){
  z <- tolower(trimws(as.character(x)))
  ifelse(z %in% c("control","ctrl") | grepl("con", z), "control", "case")
}

normalize_breed <- function(x){
  z <- tolower(trimws(as.character(x)))
  ifelse(grepl("arab|ara", z), "arabian",
         ifelse(grepl("thor|thorough", z), "thoroughbred", z))
}

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

make_dir <- function(path){
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

save_gg <- function(p, filename, w=7, h=5, dpi=300){
  ggsave(filename, plot=p, width=w, height=h, dpi=dpi)
}

# volcano plot
volcano_plot <- function(res_df,
                         title,
                         lfc_cutoff = 1,
                         fdr_cutoff = 0.05,
                         x_cap = 5,
                         y_cap = 50){
  
  df <- res_df %>%
    mutate(
      negLog10Padj = -log10(padj),
      Status = case_when(
        !is.na(padj) & padj < fdr_cutoff & log2FoldChange >  lfc_cutoff ~ "Up",
        !is.na(padj) & padj < fdr_cutoff & log2FoldChange < -lfc_cutoff ~ "Down",
        TRUE ~ "Not significant"
      ),
      # Which points exceed caps?
      x_oob = is.finite(log2FoldChange) & (log2FoldChange < -x_cap | log2FoldChange > x_cap),
      y_oob = is.finite(negLog10Padj) & (negLog10Padj > y_cap),
      any_oob = x_oob | y_oob,
      
      # "Cap" coordinates (keep points but clamp to edges)
      x_plot = pmax(pmin(log2FoldChange, x_cap), -x_cap),
      y_plot = pmin(negLog10Padj, y_cap)
    )
  
  n_x <- sum(df$x_oob, na.rm = TRUE)
  n_y <- sum(df$y_oob, na.rm = TRUE)
  n_any <- sum(df$any_oob, na.rm = TRUE)
  
  subtitle_txt <- paste0(
    "Caps applied: |log2FC| ≤ ", x_cap,
    " and -log10(FDR) ≤ ", y_cap,
    " (points beyond caps are shown on plot borders; total capped: ", n_any,
    "; x-capped: ", n_x, "; y-capped: ", n_y, ")"
  )
  
  ggplot(df, aes(x = x_plot, y = y_plot)) +
    
    # Base points (all points shown; oob points remain but are clamped)
    geom_point(aes(color = Status), size = 1.6, alpha = 0.75) +
    
    # Highlight capped points with a different shape + black outline
    geom_point(
      data = df %>% filter(any_oob),
      aes(color = Status),
      shape = 21, fill = "white", stroke = 0.5, size = 2.0, alpha = 1
    ) +
    
    scale_color_manual(
      values = c(
        "Up" = "#D73027",
        "Down" = "#4575B4",
        "Not significant" = "grey70"
      )
    ) +
    
    # Threshold lines
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff),
               linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = -log10(fdr_cutoff),
               linetype = "dashed", linewidth = 0.4) +
    
 
    scale_x_continuous(limits = c(-x_cap, x_cap)) +
    scale_y_continuous(limits = c(0, y_cap)) +
    
    labs(
      title = title,
      subtitle = subtitle_txt,
      x = paste0("log2 fold change (capped at ±", x_cap, ")"),
      y = paste0("-log10(FDR) (capped at ", y_cap, ")"),
      color = "DEG status"
    ) +
    
    theme_classic(base_size = 12) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 9)
    )
}


# PCA plot 
pca_plot <- function(vst_mat, coldata, title){
  mat <- t(vst_mat)  # samples x genes
  pcs <- prcomp(mat, scale. = FALSE)
  pvar <- (pcs$sdev^2) / sum(pcs$sdev^2)
  
  df <- data.frame(
    Sample = rownames(coldata),
    PC1 = pcs$x[,1],
    PC2 = pcs$x[,2],
    Condition = coldata$Condition
  )
  
  ggplot(df, aes(x=PC1, y=PC2, shape=Condition, color=Condition)) +
    geom_point(size=3) +
    labs(
      title=title,
      x=paste0("PC1 (", round(100*pvar[1], 1), "%)"),
      y=paste0("PC2 (", round(100*pvar[2], 1), "%)")
    ) +
    theme_classic()
}

# Heatmap: top genes; 
heatmap_top30 <- function(vst_mat, coldata, top_genes, out_png, title){
  m <- vst_mat[top_genes, , drop=FALSE]
  
  # z-score per gene
  m_z <- t(scale(t(m)))
  m_z[is.na(m_z)] <- 0
  
  
  ord <- c(which(coldata$Condition == "case"),
           which(coldata$Condition == "control"))
  m_z <- m_z[, ord, drop=FALSE]
  coldata2 <- coldata[ord, , drop=FALSE]
  
  ann <- data.frame(Condition = coldata2$Condition)
  rownames(ann) <- rownames(coldata2)
  
  png(out_png, width=1800, height=1200, res=200)
  pheatmap(
    m_z,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_col = ann,
    show_colnames = TRUE,
    show_rownames = TRUE,
    fontsize_row = 7,
    main = title
  )
  dev.off()
}

choose_design <- function(meta_b){
  tab <- table(meta_b$Horse_ID, meta_b$Condition)
  has_pairs <- any(rowSums(tab > 0) == 2)
  if (has_pairs) as.formula("~ Horse_ID + Condition") else as.formula("~ Condition")
}
