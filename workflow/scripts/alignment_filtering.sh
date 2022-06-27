#!/usr/bin/env bash

input=$1
quality=$2
threads=$3
output=$4

# Collect only deduped, mapped, paired reads of >q20
samtools idxstats ${input} | \
    cut -f 1 | \
    grep -vE 'chrM|_random|chrU|chrEBV|\*' | \
    xargs samtools view -@ $threads -f 1 -F 1284 -q $quality -o ${output} ${input}

samtools index ${output}
