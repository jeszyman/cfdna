# [[file:~/repos/cfdna-wgs/cfdna-wgs.org::*Count%20fragments%20intersecting%20windows][Count fragments intersecting windows:2]]
#!/usr/bin/env bash
input=$1
keep_bed=$2
output=$3

bedtools intersect -c \
             -a $keep_bed \
             -b $input > $output
# Count fragments intersecting windows:2 ends here
