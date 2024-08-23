#!/usr/bin/env python3
"""
Script to sample the end motifs at the 5' end of reads in a BAM file.

This script uses samtools and bedtools to extract motifs of a specified length at the 5' end
of reads in a BAM file, and then counts the occurrence of each unique motif.

Example usage:
    python3 sample_motifs.py --bam_file path/to/file.bam --output_file path/to/output.txt --motif_length 4 --threads 4 --num_reads 100000 --reference_genome path/to/reference.fasta
"""

import argparse
import itertools
import os
import logging
import sh
from collections import Counter

# --- Load Inputs --- #
# ------------------- #

def load_inputs():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam_file", type=str, required=True, help="Path to BAM file")
    parser.add_argument("--output_file", type=str, required=True, help="Path to output file")
    parser.add_argument("--motif_length", type=int, required=False, default=4, help="Length of the motif to extract (default: 4)")
    parser.add_argument("--num_reads", type=int, required=False, help="Number of reads to process (default: all)")
    parser.add_argument("--threads", type=int, required=False, default=1, help="Number of threads to use with samtools (default: 1)")
    parser.add_argument("--reference_genome", type=str, required=True, help="Path to the reference genome")
    parser.add_argument("--seed", type=int, required=False, default=12345, help="Seed for subsampling (default: 12345)")

    args = parser.parse_args()

    if not os.path.exists(args.bam_file):
        raise FileNotFoundError(f"Input BAM file not found: {args.bam_file}")

    if args.reference_genome and not os.path.exists(args.reference_genome):
        raise FileNotFoundError(f"Reference genome not found: {args.reference_genome}")

    return args

# --- Main Function --- #
# --------------------- #

def main():
    args = load_inputs()
    logging.basicConfig(level=logging.INFO)
    logging.info(f"Processing BAM file: {args.bam_file}")

    factor = calculate_sampling_factor(args.bam_file, args.num_reads)
    logging.info(f"Calculated sampling factor: {factor}")

    # Initialize motifs with zero counts for all possible motifs
    all_possible_motifs = generate_all_possible_motifs(args.motif_length)
    all_motifs = Counter({motif: 0 for motif in all_possible_motifs})

    # Apply end motif sampler to both forward and reverse flags
    forward_motifs = end_motif_sampler(args.bam_file, args.motif_length, args.num_reads, args.threads, args.reference_genome, factor, args.seed, 65)
    reverse_motifs = end_motif_sampler(args.bam_file, args.motif_length, args.num_reads, args.threads, args.reference_genome, factor, args.seed, 129)
    all_motifs.update(forward_motifs)
    all_motifs.update(reverse_motifs)

    logging.info(f"Writing results to: {args.output_file}")
    write_output(args.output_file, all_motifs)

# --- Functions --- #
# ----------------- #

def generate_all_possible_motifs(motif_length):
    """Generate all possible motifs of a given length."""
    bases = 'ACGT'
    return [''.join(p) for p in itertools.product(bases, repeat=motif_length)]

def calculate_sampling_factor(bam_file, num_reads):
    """Calculate the samtools sampling factor based on the intended number of reads."""
    f_reads = 3 * num_reads

    # Run samtools idxstats and calculate total reads
    total_reads_cmd = f"samtools idxstats {bam_file} | cut -f3 | awk '{{total += $1}} END {{print total}}'"
    try:
        total_reads = int(sh.bash('-c', total_reads_cmd).strip())
    except Exception as e:
        logging.error(f"Error calculating total reads: {e}")
        raise

    # Calculate the sampling factor
    factor = f_reads / total_reads
    return factor

def end_motif_sampler(bam_file, motif_length, num_reads, threads, reference_genome, factor, seed, flags):
    """Extract motifs from the 5' end of reads in the BAM file using the specified flags."""
    logging.info(f"Extracting motifs with motif length {motif_length} and flags {flags}.")

    cmd = f"""
    samtools view \
        --with-header \
        --min-MQ 60 \
        --require-flags {flags} \
        --subsample {factor} \
        --subsample-seed {seed} \
        --threads {threads} {bam_file} |
    bedtools bamtobed -i stdin | head -n {num_reads} |
    bedtools getfasta -bed stdin -fi {reference_genome} |
    sed "1d; n; d" | sed -E "s/(.{{{motif_length}}}).*/\\1/"
    """

    try:
        result = sh.bash('-c', cmd)
        motifs = count_motifs(result, motif_length)
    except Exception as e:
        logging.error(f"Error extracting motifs with flags {flags}: {e}")
        raise

    return motifs

def count_motifs(result, motif_length):
    """Count motifs from the result stream."""
    motifs = Counter()
    for line in result.splitlines():
        motif = line.strip()
        if len(motif) == motif_length:
            motifs[motif] += 1
    return motifs

def write_output(output_file, motifs):
    with open(output_file, 'w') as f:
        for motif, count in sorted(motifs.items()):
            f.write(f"{motif}\t{count}\n")

# --- Main Guard --- #
# ------------------ #

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
