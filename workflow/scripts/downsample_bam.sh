## Calculate the sampling factor based on the intended number of reads:
cp $1 $3

## Calculate the sampling factor based on the intended number of reads:
FACTOR=$(samtools idxstats $1 | cut -f3 | awk -v COUNT=$2 'BEGIN {total=0} {total += $1} END {print COUNT/total}')

if [[ $FACTOR > 1 ]]; then 
    echo "DS reads exceeds total for $1"
    cp $1 $3
else
    sambamba view -s $FACTOR -f bam -l 5 $1 > $3
fi
