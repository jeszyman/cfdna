#!/usr/bin/env bash

# Snakemake variables
input_bam="$1"
params_fasta="$2"
threads="${3}"
output_frag_bed="$4"

# Function
bam_to_frag(){
    # Ensure name-sorted bam file
    samtools sort -@ $threads -n -o - $1 |
    samtools fixmate -@ $threads -m -r - - |
    # Make bedpe
    bedtools bamtobed -bedpe -i - |
    # Filter any potential non-standard alignments
    awk '$1==$4 {print $0}' | awk '$2 < $6 {print $0}' |
    # Create full-fragment bed file
    awk -v OFS='\t' '{print $1,$2,$6}' |
    # Annotate with GC content and fragment length
    bedtools nuc -fi $2 -bed stdin |
    # Convert back to standard bed with additional columns
    awk -v OFS='\t' '{print $1,$2,$3,$5,$12}' |
    sed '1d' > $3
}

# Run command
bam_to_frag $input_bam \
            $params_fasta \
            $output_frag_bed
