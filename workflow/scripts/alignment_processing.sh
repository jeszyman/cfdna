#!/usr/bin/env bash


#  Note: This script is tangled from code blocks in the Emacs Org-mode file at
#  https://github.com/jeszyman/cfdna-wgs/blob/master/cfdna-wgs.org. Changes
#  made directly to this file will be overwritten upon tangle from Org-mode.


input=$1
threads=$2
output_bam=$3
output_dedup=$4
output_sort=$5
output_index=$6

sambamba view -t $threads -S -f bam $input > $output_bam
sambamba markdup -r -t $threads $output_bam $output_dedup
sambamba sort -t $threads $output_dedup -o $output_sort
sambamba index -t $threads $output_sort
