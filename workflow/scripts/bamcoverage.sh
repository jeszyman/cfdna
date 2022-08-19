#!/usr/bin/env bash

in_bam=$1
bin=$2
blacklist=$3
threads=$4
out_bg=$5

bamCoverage \
    --bam $in_bam \
    --binSize $bin \
    --blackListFileName $blacklist \
    --effectiveGenomeSize 2913022398 \
    --extendReads \
    --ignoreDuplicates \
    --ignoreForNormalization chrX \
    --normalizeUsing RPGC \
    --numberOfProcessors $threads \
    --outFileFormat bedgraph \
    --outFileName $out_bg
