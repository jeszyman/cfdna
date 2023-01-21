#!/usr/bin/env Rscript

# For unit testing
frags_tsv = "test/analysis/cfdna_wgs/frag/frag_counts.tsv"
ratios_tsv = "/home/jeszyman/mpnst/analysis/cfdna-wgs/frag/ratios.tsv"

args = commandArgs(trailingOnly = TRUE)
frags_tsv = args[1]
ratios_tsv = args[2]

# Load necessary packages
library(tidyverse)

# Load aggregate frag tsv
frags = read_tsv(frags_tsv)

# From per-position, per library short and long fragment counts, zero-centered fragment ratio
# See https://github.com/cancer-genomics/reproduce_lucas_wflow/blob/master/analysis/fig2a.Rmd

ratios =
  frags %>%
  mutate_at(vars(start, end, count), as.numeric) %>%
  # Put lib-bin short and long values on same row in order to make per-row ratios
  pivot_wider(names_from = len_class, values_from = count, values_fn = function(x) mean(x)) %>%
  mutate(fract = short/long) %>%
  select(library, chr, start, end, fract) %>%
  # Zero center by library
  group_by(library) %>%
  mutate(ratio.centered = scale(fract, scale=F)[,1])

write_tsv(ratios, file = ratios_tsv)
