#!/usr/bin/env bash
input=$1
genome_fasta=$2
output=$3

picard CollectWgsMetrics \
       INPUT=$input \
       OUTPUT=$output \
       READ_LENGTH=150 \
       REFERENCE_SEQUENCE=$genome_fasta
