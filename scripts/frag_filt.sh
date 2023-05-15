#!/usr/bin/env bash

# Steps
## Filter by absolute value of TLEN for each read
sambamba view -t $5 $1 | awk -F'\t' -v upper="$4" 'sqrt($9*$9) < upper {print $0}' | awk -F'\t' -v lower="$3" 'sqrt($9*$9) > lower {print $0}'> $2

## Restore header
sambamba view -H $1 > $6

cat $6 $2 | sambamba view -t 4 -S -f bam /dev/stdin | sambamba sort -t 4 -o $7 /dev/stdin
