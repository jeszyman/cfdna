#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
distro_dir = args[1]
healthy_libs_str = args[2]
healthy_med_file = args[3]

library(tidyverse)

healthy_libs_distros = unlist(strsplit(healthy_libs_str, " "))

read_in_gc = function(gc_csv){
  read.csv(gc_csv, header = T)
}

healthy_list = lapply(healthy_libs_distros, read_in_gc)

# Bind
healthy_all = do.call(rbind, healthy_list)

# Summarize
healthy_med =
  healthy_all %>%
  group_by(gc_strata) %>%
  summarise(med_frag_fract = median(fract_frags))

# Save
saveRDS(healthy_med, file = healthy_med_file)
