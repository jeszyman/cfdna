## Calculate the sampling factor based on the intended number of reads:
SUM=$(samtools idxstats $1 | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {print total}')
FACTOR=$(awk -v MILREADS=$2 -v SUM=$SUM 'BEGIN {print MILREADS*1000000/SUM}' )
echo "Downsampling: $1"
echo "Total reads: $SUM, Downsample target: $2 million reads, Downsample factor: $FACTOR"

if [[ $FACTOR > 1 ]]; then
    echo "DS reads exceeds total for $1"
    cp $1 $3
else
    sambamba view -s $FACTOR -f bam -l 5 -t $4 $1 > $3
    samtools index -@ $4 $3
fi
