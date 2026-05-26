# ----------------- PLOTS -----------------

volcano_plot_capped <- function(res_df,
                                title,
                                lfc_cutoff = 1,
                                fdr_cutoff = 0.05,
                                x_cap = 5,
                                y_cap = 20){
  
  df <- res_df %>%
    mutate(
      negLog10Padj = -log10(padj),
      Status = case_when(
        !is.na(padj) & padj < fdr_cutoff & log2FoldChange >  lfc_cutoff ~ "Up",
        !is.na(padj) & padj < fdr_cutoff & log2FoldChange < -lfc_cutoff ~ "Down",
        TRUE ~ "Not significant"
      ),
      x_oob = is.finite(log2FoldChange) & (log2FoldChange < -x_cap | log2FoldChange > x_cap),
      y_oob = is.finite(negLog10Padj) & (negLog10Padj > y_cap),
      any_oob = x_oob | y_oob,
      x_plot = pmax(pmin(log2FoldChange, x_cap), -x_cap),
      y_plot = pmin(negLog10Padj, y_cap)
    )
  
  n_x <- sum(df$x_oob, na.rm=TRUE)
  n_y <- sum(df$y_oob, na.rm=TRUE)
  n_any <- sum(df$any_oob, na.rm=TRUE)
  
  subtitle_txt <- paste0(
    "Caps: |log2FC|≤", x_cap, ", -log10(FDR)≤", y_cap,
    " (capped points at borders; total capped: ", n_any,
    "; x-capped: ", n_x, "; y-capped: ", n_y, ")"
  )
  
  ggplot(df, aes(x=x_plot, y=y_plot)) +
    geom_point(aes(color=Status), size=1.5, alpha=0.75) +
    geom_point(data=df %>% filter(any_oob),
               aes(color=Status),
               shape=21, fill="white", stroke=0.5, size=2.0, alpha=1) +
    scale_color_manual(values=c("Up"="#D73027","Down"="#4575B4","Not significant"="grey70")) +
    geom_vline(xintercept=c(-lfc_cutoff, lfc_cutoff), linetype="dashed", linewidth=0.4) +
    geom_hline(yintercept=-log10(fdr_cutoff), linetype="dashed", linewidth=0.4) +
    scale_x_continuous(limits=c(-x_cap, x_cap)) +
    scale_y_continuous(limits=c(0, y_cap)) +
    labs(title=title, subtitle=subtitle_txt,
         x=paste0("log2FC (capped ±", x_cap, ")"),
         y=paste0("-log10(FDR) (capped ", y_cap, ")"),
         color="DEG status") +
    theme_classic(base_size=12) +
    theme(plot.title=element_text(face="bold"),
          plot.subtitle=element_text(size=9))
}

pca_plot <- function(vst_mat, coldata, group_col, title){
  mat <- t(vst_mat)  # samples x genes
  pcs <- prcomp(mat, scale.=FALSE)
  pvar <- (pcs$sdev^2) / sum(pcs$sdev^2)
  
  df <- data.frame(
    Sample = rownames(coldata),
    PC1 = pcs$x[,1],
    PC2 = pcs$x[,2],
    Group = coldata[[group_col]]
  )
  
  ggplot(df, aes(PC1, PC2, color=Group, shape=Group)) +
    geom_point(size=3) +
    labs(title=title,
         x=paste0("PC1 (", round(100*pvar[1],1), "%)"),
         y=paste0("PC2 (", round(100*pvar[2],1), "%)"),
         color=group_col, shape=group_col) +
    theme_classic()
}

heatmap_top30 <- function(vst_mat, coldata, top_genes, out_png, title, group_col, left_level, right_level){
  m <- vst_mat[top_genes, , drop=FALSE]
  m_z <- t(scale(t(m)))
  m_z[is.na(m_z)] <- 0
  
  # fixed sample order: left group then right group
  ord <- c(which(coldata[[group_col]] == left_level),
           which(coldata[[group_col]] == right_level))
  m_z <- m_z[, ord, drop=FALSE]
  cd2 <- coldata[ord, , drop=FALSE]
  
  ann <- cd2[, intersect(c(group_col, "Breed", "Performance", "Condition"), colnames(cd2)), drop=FALSE]
  
  png(out_png, width=2000, height=1200, res=200)
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
