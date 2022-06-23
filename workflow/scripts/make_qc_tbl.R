args = commandArgs(trailingOnly = TRUE)
samstats = args[1]
read_qc_tbl = args[2]

library(dplyr)

read_qc = as_tibble(read.table(samstats, header = TRUE, sep = '\t')) %>%
  filter(grepl("dedup", Sample)) %>%
  mutate(library_id = substr(Sample,1,6)) %>%
  mutate(dedup_reads_properly_paired = reads_properly_paired) %>%
  select(library_id, dedup_reads_properly_paired)

write.table(read_qc, file = read_qc_tbl, row.names=F, sep = '\t', quote = F)
