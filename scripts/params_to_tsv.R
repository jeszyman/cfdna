#!/usr/bin/env Rscript
# ichor_params_summary_all.R
#
# Parse ALL ichorCNA *.params.txt files found recursively under --ichor-root.
# No subdir filtering. No library list filtering. Just everything.
#
# Example:
#   Rscript ichor_params_summary_all.R \
#     --ichor-root /mnt/gcs/jeszyman/projects/nf1/cfdna-cna/ichor \
#     --out summary.tsv

suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(fs)
  library(stringr)
  library(readr)
})

# ---------- CLI ----------
parser <- ArgumentParser(description = "Summarize all ichorCNA *.params.txt under a root path")
parser$add_argument("--ichor-root", type = "character", default = ".",
                    help = "Root directory to scan recursively [default: %(default)s]")
parser$add_argument("--out", type = "character", default = "-",
                    help = "Output path for TSV (use '-' for stdout)")
opt <- parser$parse_args()

# ---------- Collect all params files ----------
param_files <- dir_ls(
  path    = path.expand(opt$ichor_root),
  recurse = TRUE,
  type    = "file",
  regexp  = "\\.params\\.txt$"
)

if (length(param_files) == 0) {
  warning("No *.params.txt files found under --ichor-root")
}

# ---------- Parsers ----------
#' @title Parse a single ichorCNA params file
#' @description
#' Extracts:
#' - top-line metrics: Tumor Fraction, Ploidy
#' - GC–Map correction MAD (robust to ASCII and Unicode hyphens)
#' - fields encoded in Sample: library, frag_range, downsample, ichor_set
#' - best init (row with max loglik) if an 'init' block exists
#' @param f Path to a *.params.txt file
#' @return One-row tibble
parse_params <- function(f) {
  lines <- readr::read_lines(f)

  # Two-line TSV block header
  hdr_i <- which(stringr::str_detect(lines, "^Sample\\tTumor Fraction\\tPloidy"))[1]
  if (is.na(hdr_i) || length(lines) < hdr_i + 1) {
    return(tibble(
      library        = NA_character_,
      frag_range     = NA_character_,
      downsample     = NA_character_,
      ichor_set      = NA_character_,
      tumor_fraction = NA_real_,
      ploidy         = NA_real_,
      gc_mad         = NA_real_,
      best_init      = NA_character_,
      best_loglik    = NA_real_,
      params_file    = f
    ))
  }

  tsv <- suppressMessages(readr::read_tsv(
    paste(lines[hdr_i:(hdr_i + 1)], collapse = "\n"),
    show_col_types = FALSE
  ))

  # GC–Map correction MAD (ASCII '-', en dash '–', em dash '—')
  m <- stringr::str_match(lines, "GC[\\-–—]Map correction MAD:\\s*([0-9\\.eE\\+\\-]+)")
  gc_mad <- suppressWarnings(as.numeric(stats::na.omit(m[, 2])[1]))
  if (is.na(gc_mad)) gc_mad <- NA_real_

  # Sample-derived fields
  sample <- tsv$Sample[1]
  parts  <- strsplit(sample, "\\.")[[1]]
  if (length(parts) < 4) {
    lib <- frag_range <- downsample <- ichor_set <- NA_character_
  } else {
    lib        <- parts[1]
    frag_range <- sub("^frag_", "", parts[2])
    downsample <- sub("^ds",    "", parts[3])
    downsample <- ifelse(downsample %in% c("none","_none","NONE","_NONE","None"), "none", downsample)
    ichor_set  <- sub("^ichor_", "", parts[4])
  }

  # Optional init table → pick max loglik
  init_hdr_i <- which(stringr::str_detect(lines, "^init\\t"))[1]
  best_init   <- NA_character_
  best_loglik <- NA_real_

  if (!is.na(init_hdr_i)) {
    after <- lines[(init_hdr_i + 1):length(lines)]
    next_blank_rel <- which(stringr::str_detect(after, "^\\s*$"))[1]
    end_i <- if (is.na(next_blank_rel)) length(lines) else (init_hdr_i + next_blank_rel - 1)

    init_block <- paste(lines[init_hdr_i:end_i], collapse = "\n")
    init_tbl <- try(readr::read_tsv(init_block, show_col_types = FALSE), silent = TRUE)

    if (!inherits(init_tbl, "try-error") && nrow(init_tbl) > 0 && "loglik" %in% names(init_tbl)) {
      init_tbl$loglik <- suppressWarnings(as.numeric(init_tbl$loglik))
      init_tbl <- init_tbl[!is.na(init_tbl$loglik), , drop = FALSE]
      if (nrow(init_tbl) > 0) {
        best_row    <- init_tbl[which.max(init_tbl$loglik), , drop = FALSE]
        best_init   <- as.character(best_row$init[1])
        best_loglik <- as.numeric(best_row$loglik[1])
      }
    }
  }

  tibble::tibble(
    library        = lib,
    frag_range     = frag_range,
    downsample     = downsample,
    ichor_set      = ichor_set,
    tumor_fraction = suppressWarnings(as.numeric(tsv$`Tumor Fraction`[1])),
    ploidy         = suppressWarnings(as.numeric(tsv$Ploidy[1])),
    gc_mad         = gc_mad,
    best_init      = best_init,
    best_loglik    = best_loglik,
    params_file    = f
  )
}

#' @title Summarize a set of params files
#' @param files Character vector of params file paths
#' @return Tibble summary table, sorted
summarize_params <- function(files) {
  purrr::map_dfr(files, parse_params) %>%
    arrange(library, frag_range, downsample, ichor_set)
}

# ---------- Parse & output ----------
summary_tbl <- summarize_params(param_files)

if (opt$out == "-" || opt$out == "") {
  write_tsv(summary_tbl, file = stdout())
} else {
  write_tsv(summary_tbl, path.expand(opt$out))
}
