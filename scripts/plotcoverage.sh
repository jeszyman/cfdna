#!/usr/bin/env bash
in_bam_string=$1
blacklist=$2
threads=$3
out_raw=$4
out_plot=$5

plotCoverage \
    --bamfiles $in_bam_string \
    --blackListFileName $blacklist \
    --extendReads \
    --numberOfProcessors $threads \
    --outRawCounts $out_raw \
    --plotFile $out_plot \
    --plotFileFormat pdf \
    --skipZeros
