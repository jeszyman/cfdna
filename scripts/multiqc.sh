# [[file:~/repos/cfdna-wgs/cfdna-wgs.org::*MultiQC][MultiQC:2]]
#!/usr/bin/env bash
set -o errexit   # abort on nonzero exitstatus
set -o nounset   # abort on unbound variable
set -o pipefail  # don't hide errors within pipes

# Script variables

   input="${1}"
out_name="${2}"
 out_dir="${3}"

# Functions

multiqc $input \
        --force \
        --outdir $out_dir \
        --filename $out_name
# MultiQC:2 ends here
