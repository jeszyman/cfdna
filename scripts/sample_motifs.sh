in_bam="${1}"
in_fasta="${2}"
n_motif="${3}"
n_reads="${4}"
seed="${5}"
threads="${6}"
out_merged="${7}"

main(){
    forward_motif \
        $in_bam \
        $seed \
        $threads \
        $n_reads \
        $in_fasta \
        $n_motif > $out_merged
    reverse_motif \
        $in_bam \
        $seed \
        $threads \
        $n_reads \
        $in_fasta \
        $n_motif >> $out_merged
}

#########1#########2#########3#########4#########5#########6#########7#########8
forward_motif(){
    #
    local in_bam="${1}"
    local seed="${2}"
    local threads="${3}"
    local n_reads="${4}"
    local in_fasta="${5}"
    local n_motif="${6}"
    #
    # Calculate the samtools sampling factor based on the intended number of
    # reads. This will be 2x the n_read input plus some margin of error as
    # the next set will only to forward reads.
    f_reads=$(( 3*$n_reads ))
    factor=$(samtools idxstats $in_bam |
                 cut -f3 |
                 awk -v nreads=$f_reads 'BEGIN {total=0} {total += $1} END {print nreads/total}')
    #
    # Take first read in mapped, paired, with normal FS orientation.
    # View perfect matching reads (for BWA), first in pair.
    samtools view \
             --with-header \
             --min-MQ 60 \
             --require-flags 65 \
             --subsample $factor \
             --subsample-seed $seed \
             --threads $threads $in_bam |
        # Fetch reference for n reads
        bedtools bamtobed -i stdin | head -n $n_reads |
        bedtools getfasta -bed stdin -fi $in_fasta |
        # Sed magic to extract motifs from fasta
        sed "1d; n; d" | sed -E "s/(.{$n_motif}).*/\1/"
}

#########1#########2#########3#########4#########5#########6#########7#########8
reverse_motif(){
    #
    local in_bam="${1}"
    local seed="${2}"
    local threads="${3}"
    local n_reads="${4}"
    local in_fasta="${5}"
    local n_motif="${6}"
    #
    # Calculate the samtools sampling factor based on the intended number of
    # reads. This will be 2x the n_read input plus some margin of error as
    # the next set will only to forward reads.
    f_reads=$(( 3*$n_reads ))
    factor=$(samtools idxstats $in_bam |
                 cut -f3 |
                 awk -v nreads=$f_reads 'BEGIN {total=0} {total += $1} END {print nreads/total}')
    #
    # Take SECOND read in mapped, paired, with normal FS orientation.
    # View perfect matching reads (for BWA).
    samtools view \
             --with-header \
             --min-MQ 60 \
             --require-flags 129 \
             --subsample $factor \
             --subsample-seed $seed \
             --threads $threads $in_bam |
        # Fetch reference for n reads
        bedtools bamtobed -i stdin | head -n $n_reads |
        bedtools getfasta -bed stdin -fi $in_fasta |
        # Sed magic to extract motifs from fasta
        sed "1d; n; d" | sed -E "s/.*(.{$n_motif})/\1/" |
        # Generate reverse compliment
        tr ACGT TGCA | rev
}

main "$@"
