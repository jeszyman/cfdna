#!/usr/bin/env bash
in_fasta="${1}"
out_prefix="${2}"

bwa index -p $out_prefix $in_fasta
