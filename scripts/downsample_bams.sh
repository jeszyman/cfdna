#!/usr/bin/env bash
# For unit testing
# in_bam=test/analysis/cfdna_wgs/bams/lib001_filt.bam
# milreads=0.001
# checker=test/tmp/lib001_ds0.001.txt
# outdir=test/analysis/cfdna_wgs/bams
# suffix=_filt.bam
# threads=4

in_bam="${1}"
milreads="${2}"
outdir="${3}"
suffix="${4}"
threads="${5}"

downsample(){
    # Derived variables
    milreads_full=$(awk -v milreads="${milreads}" 'BEGIN{milreads_full=(1000000*milreads); print milreads_full}')
    factor=$(samtools idxstats $in_bam |
                 cut -f3 |
                 awk -v count=$milreads_full 'BEGIN {total=0} {total += $1} END {print count/total}')
    base=$(basename -s $suffix $in_bam)
    out_bam=${outdir}/${base}_ds${milreads}.bam
    #
    # Downsample
    if [[ $factor < 1 ]]; then
    samtools view -s $factor -b -@ $threads $in_bam > $out_bam
    fi
}

downsample $in_bam $milreads $suffix
