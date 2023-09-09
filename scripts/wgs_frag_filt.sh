#!/usr/bin/env bash

#########################################
###   Filter Bam By Fragment Length   ###
#########################################

inbam="${1}"
nohead="${2}"
min="${3}"
max="${4}"
threads="${5}"
onlyhead="${6}"
outbam="${7}"

# Filter by absolute value of TLEN for each read
samtools view -@ $threads $inbam |
    awk -F'\t' -v upper="$max" 'sqrt($9*$9) < upper {print $0}' |
    awk -F'\t' -v lower="$min" 'sqrt($9*$9) > lower {print $0}'> $nohead

# Restore header
samtools view -@ $threads --header-only $inbam > $onlyhead


cat $onlyhead $nohead |
    samtools view -@ $threads --bam /dev/stdin |
    samtools sort -@ $threads -o $outbam /dev/stdin
