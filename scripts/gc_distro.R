#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
bed_file = args[1]
distro_file = args[2]

library(tidyverse)

# Read in modified bed
bed = read.table(bed_file, sep = '\t')
names(bed) = c("chr","start","end","gc_raw","len")

# Generate distribution csv
distro =
  bed %>%
  # Round GC
  mutate(gc_strata = round(gc_raw, 2)) %>%
  # Count frags per strata
  count(gc_strata) %>%
  # Get fraction frags
  mutate(fract_frags = n/sum(n)) %>% mutate(library_id = gsub("_frag.bed", "", gsub("^.*lib", "lib", bed_file))) %>%
  select(library_id,gc_strata,fract_frags) %>%
  write.csv(file = distro_file, row.names = F)
