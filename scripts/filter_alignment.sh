#!/usr/bin/env bash

input=$1
threads=$2
output=$3

# Filter to reads that are
#  - Excluding any unmapped, not primary alignment, or duplicates
#  - Only MAPQ > 20
# DO NOT restrict to "proper pairs"- this clips long cfDNA fragments!

samtools view -@ $threads -b -F 1284 -h -q 20 -o $output $input

samtools index ${output}
