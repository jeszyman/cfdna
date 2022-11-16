#!/usr/bin/env bash

# For unit test
# singularity shell ~/sing_containers/cfdna_wgs.1.0.0.sif
# blacklist="test/inputs/hg38-blacklist.v2.bed.gz"
# chrom_sizes="test/inputs/hg38.chrom.sizes"
# auto_bed="/tmp/test.bed"
# keep_bed="/tmp/keep.bed"

blacklist="${1}"
chrom_sizes="${2}"
auto_bed="${3}"
keep_bed="${4}"

# Make autosome bed from chrom_sizes
cat $chrom_sizes | grep -v _ | grep chr[0-9] | awk -v OFS='\t' '{ print $1, 0, $2}' > $auto_bed

# Filter autosome bed by blacklist
bedtools subtract -a $auto_bed -b $blacklist > $keep_bed
