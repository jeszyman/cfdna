#! /bin/bash

input_read1=$1
input_read2=$2
adapter_fasta=$3
threads=$4
output_read1=$5
output_read1_unpr=$6
output_read2=$7
output_read2_unpr=$8
report_name=$9

module load fastp # 0.23.2 on biowulf

fastp --in1 $input_read1 \
      --in2 $input_read2 \
      --out1 $output_read1 \
      --unpaired1 $output_read1_unpr \
      --out2 $output_read2 \
      --unpaired2 $output_read2_unpr \
      --thread $threads \
      --json $report_name.json \
      --html $report_name.html \
      --length_required 20 \
      --cut_front --cut_front_mean_quality 10 \
      --cut_tail --cut_tail_mean_quality 10 \
      --adapter_fasta $adapter_fasta
