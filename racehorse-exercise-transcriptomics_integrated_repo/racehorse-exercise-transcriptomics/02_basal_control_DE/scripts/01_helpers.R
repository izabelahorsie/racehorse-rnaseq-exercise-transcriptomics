suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(pheatmap)
  library(readxl)
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

save_gg <- function(p, filename, w=7, h=5, dpi=300){
  ggsave(filename, plot=p, width=w, height=h, dpi=dpi)
}
