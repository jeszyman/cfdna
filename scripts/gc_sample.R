bash [[file:~/repos/cfdna-wgs/cfdna-wgs.org::*Sample%20fragments%20by%20healthy%20GC%20proportions][Sample fragments by healthy GC proportions:2]]
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
healthy_med = args[1]
frag_file = args[2]
sampled_file = args[3]

library(tidyverse)

healthy_fract = readRDS(healthy_med)
frag_file = read.table(frag_file, sep = '\t', header = F)

frag_bed = frag_file
names(frag_bed) = c("chr", "start", "end", "gc_raw", "len")

frag = frag_bed %>%
  # Round off the GC strata
  mutate(gc_strata = round(gc_raw, 2)) %>%
  # Join the median count of fragments per strata in healthies
  # Use this later as sampling weight
  left_join(healthy_fract, by = "gc_strata")

# Determine frags to sample by counts in strata for which healthies had highest count
stratatotake = frag$gc_strata[which.max(frag$med_frag_fract)]
fragsinmaxstrata = length(which(frag$gc_strata == stratatotake))
fragstotake = round(fragsinmaxstrata/stratatotake)

sampled = frag %>%
  filter(!is.na(med_frag_fract)) %>%
  slice_sample(., n = nrow(.), weight_by = med_frag_fract, replace = T) %>% select(chr, start, end, len, gc_strata)

write.table(sampled, sep = "\t", col.names = F, row.names = F, quote = F, file = sampled_file)
bash Sample fragments by healthy GC proportions:2 ends here
