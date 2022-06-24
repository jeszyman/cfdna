#!/usr/bin/env bash

#  Note: This script is tangled from code blocks in the Emacs Org-mode file at
#  https://github.com/jeszyman/cfdna-wgs/blob/master/cfdna-wgs.org. Changes
#  made directly to this file will be overwritten upon tangle from Org-mode.


input_read1=$1
input_read2=$2
params_adapter_fasta=$3
threads=$4
output_read1=$5
output_read1_unpr=$6
output_read2=$7
output_read2_unpr=$8
log_int=$9

trimmomatic PE \
            -threads $threads \
            -trimlog $log_int \
            $input_read1 $input_read2 \
            $output_read1 $output_read1_unpr \
            $output_read2 $output_read2_unpr \
            ILLUMINACLIP:$params_adapter_fasta:2:30:10 \
            LEADING:10 TRAILING:10 MAXINFO:50:0.97 MINLEN:20
