#!/usr/bin/env bash

# For unit testing
# in_bam="test/analysis/cfdna_wgs_bams/lib001_filt.bam"
# out_bam=/tmp/test.bam
# milreads="0.0041"

in_bam=$1
milreads="$2"
out_bam=$3

reads=$(echo |awk -v var1=$milreads '{ print 1000000*var1 }')

## Calculate the sampling factor based on the intended number of reads:

FACTOR=$(samtools idxstats $in_bam | cut -f3 | awk -v COUNT=$reads 'BEGIN {total=0} {total += $1} END {print COUNT/total}')

if [[ $FACTOR > 1 ]]; then
    echo "DS reads exceeds total for $in_bam"
else
    sambamba view -s $FACTOR -f bam -l 5 $in_bam > $out_bam
fi
