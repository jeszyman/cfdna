#!/usr/bin/env bash
set -o errexit
set -o nounset
set -o pipefail

print_usage() {
    cat <<EOF
USAGE: cfdna_frag_filt.sh <IN_BAM> <NOHEAD_OUT> <MIN_TLEN> <MAX_TLEN> <THREADS> <HEADER_OUT> <FINAL_BAM> <TMPDIR>

DESCRIPTION:
  Filters BAM by absolute TLEN (template length) and outputs a sorted, filtered BAM.

REQUIRED ARGUMENTS:
  <IN_BAM>        Input BAM file
  <NOHEAD_OUT>    Output path for filtered alignments without header
  <MIN_TLEN>      Minimum TLEN (exclusive)
  <MAX_TLEN>      Maximum TLEN (exclusive)
  <THREADS>       Number of threads for samtools
  <HEADER_OUT>    Output path for BAM header
  <FINAL_BAM>     Output path for final sorted BAM
  <TMPDIR>        Temporary directory for samtools sort

EXAMPLE:
  cfdna_frag_filt.sh sample.bam filtered.sam 100 400 8 header.sam filtered_sorted.bam /tmp

EOF
}

parse_args() {
    if [[ "$#" -ne 8 ]] || [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
        print_usage
        exit 1
    fi

    in_bam="$1"
    nohead="$2"
    min_tlen="$3"
    max_tlen="$4"
    threads="$5"
    header_out="$6"
    final_bam="$7"
    tmpdir="$8"
}

main() {
    # Filter reads by |TLEN| between min and max
    samtools view -@ "$threads" "$in_bam" |
        awk -F'\t' -v min="$min_tlen" -v max="$max_tlen" '{
            len = sqrt($9 * $9);
            if (len > min && len < max) print
        }' > "$nohead"

    # Extract header
    samtools view -@ "$threads" --header-only "$in_bam" > "$header_out"

    # Combine and sort
    cat "$header_out" "$nohead" |
        samtools view -@ "$threads" --bam /dev/stdin |
        samtools sort -@ "$threads" -T "$tmpdir/sorttemp" -o "$final_bam" /dev/stdin
}

parse_args "$@"
main
