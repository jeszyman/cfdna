#!/usr/bin/env bash
set -o errexit
set -o nounset
set -o pipefail

print_usage() {
    cat <<EOF
USAGE: cfdna_downsample.sh <IN_BAM> <MILLION_READS>

DESCRIPTION:
  Downsamples BAM to an approximate target number of reads (in millions) using samtools -s.

REQUIRED ARGUMENTS:
  <IN_BAM>          Input BAM file
  <MILLION_READS>   Target number of reads, in millions (e.g. 10 for ~10 million reads)

EXAMPLE:
  cfdna_downsample.sh sample_filt.bam 10

NOTE:
  Output BAM will be named by replacing '_filt.bam' with '_ds10.bam'
EOF
}

parse_args() {
    if [[ "$#" -ne 2 ]] || [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
        print_usage
        exit 1
    fi

    in_bam="$1"
    million_reads="$2"
}

main() {
    reads=$(awk -v m="$million_reads" 'BEGIN { print m * 1000000 }')

    out_bam=$(echo "$in_bam" | sed 's/_filt.bam/_ds10.bam/g')

    factor=$(samtools idxstats "$in_bam" | cut -f3 |
        awk -v count="$reads" 'BEGIN { total = 0 } { total += $1 } END { print count / total }')

    samtools view -@ 4 -s "$factor" -b "$in_bam" > "$out_bam"
    samtools index "$out_bam"
}

parse_args "$@"
main
