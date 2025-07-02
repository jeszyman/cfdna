#!/usr/bin/env bash
set -o errexit
set -o nounset
set -o pipefail

print_usage() {
    cat <<EOF
USAGE: bwa_index.sh <INPUT_FASTA_GZ> <OUTPUT_FASTA> <OUTPUT_FAI> <OUTPUT_BED> <BWA_PREFIX> <LOG_FILE>

DESCRIPTION:
  Prepare a BWA index from a gzipped FASTA input, including .fai indexing and BED file generation.

REQUIRED ARGUMENTS:
  <INPUT_FASTA_GZ>   Path to input .fa.gz file
  <OUTPUT_FASTA>     Ungzipped target FASTA output path
  <OUTPUT_FAI>       Output .fai file (samtools faidx)
  <OUTPUT_BED>       Output BED file for primary chromosomes
  <BWA_PREFIX>       Prefix for BWA index files
  <LOG_FILE>         Log file for BWA index command

OPTIONS:
  -h, --help         Show this help message

EXAMPLE:
  bwa_index.sh hg38.fa.gz ref/hg38/hg38.fa ref/hg38/hg38.fa.fai ref/hg38/hg38.primary.bed ref/hg38/hg38 ref/hg38/bwa_index.log
EOF
}

main() {
    parse_args "$@"
    make_index_dir
    decompress_fasta
    index_fasta
    generate_bed
    run_bwa_index
}

parse_args() {
    if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
        print_usage
        exit 0
    fi

    if [[ $# -ne 6 ]]; then
        echo "Error: Expected 6 arguments."
        print_usage
        exit 1
    fi

    declare -g input_fasta_gz="$1"
    declare -g output_fasta="$2"
    declare -g output_fai="$3"
    declare -g output_bed="$4"
    declare -g bwa_prefix="$5"
    declare -g log_file="$6"
}

make_index_dir() {
    mkdir -p "$(dirname "$output_fasta")"
}

decompress_fasta() {
    echo "Decompressing $input_fasta_gz to $output_fasta"
    gzip -dc "$input_fasta_gz" > "$output_fasta"
}

index_fasta() {
    echo "Indexing $output_fasta with samtools faidx"
    samtools faidx "$output_fasta"
}

generate_bed() {
    echo "Generating BED file: $output_bed"
    cut -f1,2 "$output_fai" \
        | grep -E '^chr([1-9]|1[0-9]|2[0-2]|X|Y)\s' \
        | awk '{print $1 "\t0\t" $2}' > "$output_bed"
}

run_bwa_index() {
    echo "Running bwa index: $bwa_prefix"
    bwa index -p "$bwa_prefix" "$output_fasta" > "$log_file" 2>&1
}

main "$@"
