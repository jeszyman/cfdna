## fastqc_input="test/qc/all_qc_data/multiqc_fastqc.txt"
## samstats_input="test/qc/all_qc_data/multiqc_samtools_stats.txt"
## flagstats_input="test/qc/all_qc_data/multiqc_samtools_flagstat.txt"
## picard_input="test/qc/all_qc_data/multiqc_picard_wgsmetrics.txt"
## deeptools_input="test/qc/all_frag.tsv"

args = commandArgs(trailingOnly = TRUE)
fastqc_input = args[1]
samstats_input = args[2]
flagstats_input = args[3]
picard_input = args[4]
deeptools_input = args[5]
readqc_out_tbl = args[6]

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
  mutate(library = substr(Sample, 1, 6))

flagstats = as_tibble(read.table(flagstats_input, header = TRUE, sep = '\t', stringsAsFactors = FALSE)) %>%
  mutate(library = substr(Sample, 1, 6))

deeptools = as_tibble(read.table(deeptools_input, header = FALSE, sep = '\t', stringsAsFactors = FALSE))
colnames(deeptools)=c("frag_len","frag_count","file")
deeptools = deeptools %>%
  mutate(library = substr(file, nchar(file) -9, nchar(file) -4)) %>%
  mutate(frag_len = sub("^", "frag_len", frag_len)) %>%
  select(library, frag_len, frag_count) %>%
  pivot_wider(
    names_from = frag_len,
    values_from = frag_count)

picard = as_tibble(read.table(picard_input, header = TRUE, sep = '\t', stringsAsFactors = FALSE)) %>%
  mutate(library = Sample)

readqc = fastqc %>%
  left_join(samstats, by = "library") %>%
  left_join(flagstats, by = "library") %>%
  left_join(deeptools, by = "library") %>%
  left_join(picard, by = "library")

write.table(readqc, file = readqc_out_tbl, row.names = F, sep = '\t', quote = F)
