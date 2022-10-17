#!/usr/bin/env bash
set -o errexit   # abort on nonzero exitstatus
set -o nounset   # abort on unbound variable
set -o pipefail  # don't hide errors within pipes

# Script variables

input="${1}"
log="${2}"
output_hist="${3}"
output_raw="${4}"
blacklist="${5}"
threads="${6}"


bamPEFragmentSize --bamfiles $input \
                  --numberOfProcessors $threads \
                  --blackListFileName $blacklist \
                  --histogram $output_hist \
                  --maxFragmentLength 1000 \
                  --outRawFragmentLengths $output_raw
