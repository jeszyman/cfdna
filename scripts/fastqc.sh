# [[file:~/repos/cfdna-wgs/cfdna-wgs.org::*FastQC][FastQC:2]]
#!/usr/bin/env bash
set -o errexit   # abort on nonzero exitstatus
set -o nounset   # abort on unbound variable
set -o pipefail  # don't hide errors within pipes

# Script variables
input="${1}"
outdir="${2}"
threads="${3}"

# Functions
fastqc  --outdir $outdir \
        --quiet \
        --threads $threads $input
# FastQC:2 ends here
