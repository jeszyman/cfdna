# [[file:~/repos/cfdna-wgs/cfdna-wgs.org::*Make%20GC%20and%20mappability%20restricted%20bins][Make GC and mappability restricted bins:2]]
gc5mb="${1}"
blklist="${2}"
keep="${3}"

bedtools intersect -a $gc5mb -b $blklist -v -wa |
    grep -v _ |
    awk '{ if ($4 >= 0.3) print $0 }' > $keep
# Make GC and mappability restricted bins:2 ends here
