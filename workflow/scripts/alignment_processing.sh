#!/usr/bin/env bash

<#bash_preamble#>

input=$1
threads=$2
output_bam=$3
output_dedup=$4
output_sort=$5
output_index=$6

sambamba view -t $threads -S -f bam $input > $output_bam
sambamba markdup -r -t $threads $output_bam $output_dedup
sambamba sort -t $threads $output_dedup -o $output_sort
sambamba index -t $threads $output_sort
