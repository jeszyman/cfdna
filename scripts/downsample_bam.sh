#!/usr/bin/env bash

in_bam=$1
milreads="$2"

reads=$(echo |awk -v var1=$milreads '{ print 1000000*var1 }')

out_bam=$(echo $in_bam | sed 's/_filt.bam/_ds10.bam/g')

## Calculate the sampling factor based on the intended number of reads:

FACTOR=$(samtools idxstats $in_bam | cut -f3 | awk -v COUNT=$reads 'BEGIN {total=0} {total += $1} END {print COUNT/total}')

samtools view -@ 4 -s $FACTOR -b $in_bam > $out_bam

samtools index $out_bam
