gc5mb="${1}"
blklist="${2}"
keep="${3}"

bedtools intersect -a $gc5mb -b $blklist -v -wa |
    grep -v _ |
    awk '{ if ($4 >= 0.3) print $0 }' > $keep
