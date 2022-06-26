#fastqc_input="test/qc/all_qc_data/multiqc_fastqc.txt"
#samstats_input="test/qc/all_qc_data/multiqc_samtools_stats.txt"
#flagstats_input="test/qc/all_qc_data/multiqc_samtools_flagstat.txt"

args = commandArgs(trailingOnly = TRUE)
fastqc_input = args[1]
samstats_input = args[2]
flagstats_input = args[3]
readqc_out_tbl = args[4]

library(tidyverse)

fastqc = as_tibble(read.table(fastqc_input, header = TRUE, sep = '\t', stringsAsFactors = FALSE)) %>%
  mutate(library = substr(Filename,1,6)) %>%
  mutate(read = ifelse(grepl("R1", Filename), "read1", "read2")) %>%
  mutate(fastq_processing = ifelse(grepl("proc", Filename), "processed", "raw")) %>%
  select(!c(Sample,File.type,Encoding)) %>%
  pivot_wider(
    names_from = c(read,fastq_processing),
    values_from = !c(library,read,fastq_processing))

samstats = as_tibble(read.table(samstats_input, header = TRUE, sep = '\t', stringsAsFactors = FALSE)) %>%
  mutate(library = substr(Sample, 1, 6)) %>%
  mutate(bam_processing = ifelse(grepl("dedup",Sample), "dedup", "raw")) %>%
  pivot_wider(
    names_from = bam_processing,
    values_from = !c(library, bam_processing))

flagstats = as_tibble(read.table(flagstats_input, header = TRUE, sep = '\t', stringsAsFactors = FALSE)) %>%
  mutate(library = substr(Sample, 1, 6)) %>%
  mutate(bam_processing = ifelse(grepl("dedup",Sample), "dedup", "raw")) %>%
  pivot_wider(
    names_from = bam_processing,
    values_from = !c(library, bam_processing))

readqc = fastqc %>% left_join(samstats, by = "library") %>% left_join(flagstats, by = "library")

write.table(readqc, file = readqc_out_tbl, row.names = F, sep = '\t', quote = F)
