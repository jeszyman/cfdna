#!/usr/bin/env bash
output=$1
declare -a array2=$2

if [ -f $output ]; then \rm $output; fi

for file in ${array2[@]}; do
    awk '{{print FILENAME (NF?"\t":"") $0}}' $file |
        sed 's/^.*lib/lib/g' |
        sed 's/_.*_/\t/g' |
        sed 's/\.bed//g' >> $output
done
