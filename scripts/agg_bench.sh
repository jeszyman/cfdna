# For unit testing
indir="test/benchmark"
output="test/analysis/qc/bench_agg.tsv"

if [ -f $output ]; then rm $output; fi

for file in $indir/*
do
    base=$(basename $file)
    cat $file | awk -v OFS='\t' -v var=$base 'NR>1 {print var,$0}' >> $output
done

sed -i '1i\process\tfloat_sec\trun_time\tmax_rss\tmax_vms\tmax_uss\tmax_pss\tio_in\tio_out\tmean_load\tcpu_time' $output
