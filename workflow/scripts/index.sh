#!/usr/bin/env bash
input=$1
out_prefix=$2
bwa index -p {out_prefix} {input}
