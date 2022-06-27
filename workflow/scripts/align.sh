#!/usr/bin/env bash

input_ref=$1
input_r1=$2
input_r2=$3
threads=$4
output_sort=$5

bwa mem -M -t $threads \
    $input_ref \
    $input_r1 \
    $input_r2 |
    samblaster -M |
    samtools view -@ $threads -Sb - -o - |
    samtools sort -@ $threads - -o $output_sort
samtools index -@ threads $output_sort
