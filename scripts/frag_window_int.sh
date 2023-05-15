#!/usr/bin/env bash
input=$1
keep_bed=$2
output=$3

bedtools intersect -c \
             -a $keep_bed \
             -b $input > $output
