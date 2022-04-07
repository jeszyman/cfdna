#!/usr/bin/env bash
trimmomatic PE \
            -threads 2 \
            ~/repos/mpnst-preprocessing/test/fastq/mpnst1_R1.fastq.gz \
            ~/repos/mpnst-preprocessing/test/fastq/mpnst1_R2.fastq.gz \
            /tmp/test_R1_paired.fastq.gz /tmp/test_R2_paired.fastq.gz \
            /tmp/test_R1_unpaied.fastq.gz /tmp/test_R2_unpaired.fastq.gz \
            ILLUMINACLIP:~/repos/mpnst-preprocessing/test/inputs/TruSeq3-PE.fa:2:30:10
