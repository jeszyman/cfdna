#!/usr/bin/env bash

input=$1
keepbed=$2
threads=$3
output=$4

# Filter to reads that are
#  - Only mapped in proper pairs (-f 3)
#  - Excluding any unmapped, not primary alignment, or duplicates
#  - Only mapped to regions in the keep.bed file (-L $bed) (autosomes not in blacklist)
#  - Only MAPQ > 20

samtools view -@ $threads -b -f 3 -F 1284 -h -L $keepbed -M -q 20 -o $output $input

samtools index ${output}
