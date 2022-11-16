#!/usr/bin/env Rscript
#fastqc_input="test/qc/all_cfdna_wgs_data/multiqc_fastqc.txt"
## samstats_input="test/qc/all_cfdna_wgs_data/multiqc_samtools_stats.txt"
## flagstats_input="test/qc/all_cfdna_wgs_data/multiqc_samtools_flagstat.txt"
## picard_input="test/qc/all_cfdna_wgs_data/multiqc_picard_wgsmetrics.txt"
deeptools_frag_input="test/qc/deeptools_frag_lengths.txt"

args = commandArgs(trailingOnly = TRUE)
fastqc_input = args[1]
samstats_input = args[2]
flagstats_input = args[3]
picard_input = args[4]
deeptools_frag_input = args[5]
deeptools_cov_input = args[6]
readqc_out_tbl = args[7]
frag_len_out_tbl = args[8]

library(tidyverse)

process_multiqc_fastqc = function(multiqc_fastqc_input){
  as_tibble(read.table(multiqc_fastqc_input, header = TRUE, sep = '\t', stringsAsFactors = FALSE)) %>%
  mutate(library = substr(Filename,1,6)) %>%
  mutate(read = ifelse(grepl("R1", Filename), "read1", "read2")) %>%
  mutate(fastq_processing = gsub("_.*$","",substr(Sample, 8, length(Sample)))) %>%
  select(!c(Sample,File.type,Encoding)) %>%
  pivot_wider(
    names_from = c(read,fastq_processing),
    values_from = !c(library,read,fastq_processing))
}

fastqc = process_multiqc_fastqc(fastqc_input)

process_multiqc_samfile = function(multiqc_samfile){
  as_tibble(read.table(multiqc_samfile, header = TRUE, sep = '\t', stringsAsFactors = FALSE)) %>%
  mutate(library = substr(Sample, 1, 6)) %>%
  mutate(bam_processing = gsub("_.*$","",substr(Sample, 8, length(Sample)))) %>%
  select(!c(Sample)) %>%
  pivot_wider(
    names_from = c(bam_processing),
    values_from = !c(library, bam_processing))
}

samstats = process_multiqc_samfile(samstats_input)
flagstats = process_multiqc_samfile(flagstats_input)

deeptools_frag = read_tsv(deeptools_frag_input, col_names = c("frag_len","frag_count","file"), skip = 1) %>%
  filter(frag_len < 500) %>%
  mutate(library = substr(gsub("^.*lib", "lib", file), 1,6)) %>%
  mutate(frag_len = sub("^", "frag_len", frag_len)) %>%
  select(library, frag_len, frag_count) %>%
  pivot_wider(
    names_from = frag_len,
    values_from = frag_count)

picard = as_tibble(read.table(picard_input, header = TRUE, sep = '\t', stringsAsFactors = FALSE)) %>%
  mutate(library = Sample)

deeptools_cov = read_tsv(deeptools_cov_input, skip = 1) %>%
  pivot_longer(!c(`#'chr'`, `'start'`,`'end'`), names_to = "file", values_to = "cnt") %>%
  rename(chr = `#'chr'`,
         start = `'start'`,
         end = `'end'`) %>%
  mutate(library = substr(file, 2, 7)) %>%
  group_by(library) %>%
  summarise(
    mean_cov = mean(cnt),
    median_cov = median(cnt),
            )

readqc = fastqc %>%
  left_join(samstats, by = "library") %>%
  left_join(flagstats, by = "library") %>%
  left_join(deeptools_frag, by = "library") %>%
  left_join(picard, by = "library") %>%
  left_join(deeptools_cov, by = "library")

write.table(readqc, file = readqc_out_tbl, row.names = F, sep = '\t', quote = F)

all_frag_len = data.frame(frag_len = 1:500)

frag_len =
  readqc %>% select(starts_with("frag_len") | matches("library")) %>%
  pivot_longer(!library, names_to = "frag_len", values_to = "count") %>%
  mutate(frag_len = as.numeric(gsub("frag_len","",frag_len))) %>%
  mutate(count = as.numeric(count)) %>%
  pivot_wider(names_from = library, values_from = count) %>%
  right_join(all_frag_len) %>% arrange(frag_len) %>%
  replace(is.na(.), 0)

write_tsv(frag_len, file = frag_len_out_tbl)
