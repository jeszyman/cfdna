# [[file:~/repos/cfdna-wgs/cfdna-wgs.org::*Remove%20PCR%20duplicates][Remove PCR duplicates:2]]
#!/usr/bin/env bash
set -o errexit   # abort on nonzero exitstatus
set -o nounset   # abort on unbound variable
set -o pipefail  # don't hide errors within pipes

# Script variables
raw_bam="${1}"
dedup_bam="${2}"
threads="${3}"

samtools sort -@ $threads -n -o - $raw_bam |
    samtools fixmate -m - - |
    samtools sort -@ $threads -o - - |
    samtools markdup -@ $threads -r - $dedup_bam
samtools index $dedup_bam
# Remove PCR duplicates:2 ends here
