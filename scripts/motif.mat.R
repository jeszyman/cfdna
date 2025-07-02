#!/usr/bin/env Rscript

########################################
###   Make End Motif Single Matrix   ###
########################################

args = commandArgs(trailingOnly = TRUE)
motif_str = args[1]
motif_tsv = args[2]

# Load required packages, data, and functions
library(tidyverse)

# Define possible 4-mer motifs
possible_motifs =
  expand.grid(rep(list(c('A', 'G', 'T', 'C')), 4)) %>%
  as_tibble() %>%
  mutate(motif = paste0(Var1,Var2,Var3,Var4)) %>%
  select(motif) %>% arrange(motif)
possible_motifs

# Define motif files list
#motif_str = "~/mpnst/analysis/frag/motifs/lib005_motifs.txt ~/mpnst/analysis/frag/motifs/lib507_motifs.txt"
(motif_files = strsplit(motif_str, " ")[[1]])
(names(motif_files) = substr(gsub("^.*lib","lib",motif_files), 1, 6))


#(motif_files = list.files(motif_samples_dir, full.names = TRUE, pattern = "^lib.*motifs.txt"))
#(names(motif_files)=substr(list.files(motif_samples_dir, pattern = "^lib.*motifs.txt"),1,6))

# Make per-libary motif frequencies
ingest_motif = function(motif_file){
  read_tsv(motif_file,
           col_names = c("motif")) %>%
    group_by(motif) %>%
    summarise(count = n()) %>%
    mutate(fract = count/sum(count)) %>%
    select(motif, fract)
}

motif_tibs = lapply(motif_files, ingest_motif)

# Make single matrix tsv
motifs = bind_rows(motif_tibs, .id = "library") %>% pivot_wider(names_from = library, values_from = fract) %>% filter(motif %in% possible_motifs$motif)

motifs %>% write_tsv(., motif_tsv)
