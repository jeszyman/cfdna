# [[file:~/repos/cfdna-wgs/cfdna-wgs.org::*deepTools%20bamCoverage][deepTools bamCoverage:2]]
#!/usr/bin/env bash

in_bam=$1
bin=$3
blacklist=$4
threads=$5
out_bg=$2

bamCoverage \
    --bam $in_bam \
    --binSize $bin \
    --blackListFileName $blacklist \
    --effectiveGenomeSize 2913022398 \
    --extendReads \
    --ignoreDuplicates \
    --ignoreForNormalization chrX \
    --normalizeUsing RPGC \
    --numberOfProcessors $threads \
    --outFileFormat bedgraph \
    --outFileName $out_bg
# deepTools bamCoverage:2 ends here
