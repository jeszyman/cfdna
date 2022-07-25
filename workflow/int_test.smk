##################################################################
###   Integration testing snakefile for WGS cfDNA Processing   ###
##################################################################

import pandas as pd
import re
container: config["container"]

# Import libraries table to pandas
libraries = pd.read_table(config["data_dir"] + "/inputs/libraries.tsv")

# Pull library ids out of libraries table
LIBRARY_IDS = list(libraries.library.unique())

# List of downsampling values in millions of reads
MILREADS = config["MILREADS"]


# Function acts on read_qc, generated in the workflow, to select libraries for
# downsampling. Notice library 2 does not downsample because it already has
# fewer than 3000 reads. Best practice for real data would be to use the
# MILREADS value in lieu of a specified number here.

def get_ds_candidates(wildcards):
    read_qc = pd.read_table(checkpoints.make_qc_tbl.get().output[0])
    test=read_qc.library[read_qc.reads_properly_paired > 3000].tolist()
    return expand(
	config["data_dir"] + "/bam/ds/{library_id}_ds{milreads}.bam",
        library_id=test, milreads = MILREADS)

# Makes the name bwa index directory from the config genome fasta
#  e.g. test/inputs/chr8.fa will make test/ref/chr8
genome_ref = config["genome_fasta"]
genome_ref = re.sub("inputs", lambda x: 'ref', genome_ref)
genome_ref = re.sub("\..*$", lambda x: '', genome_ref)

#########1#########2#########3#########4#########5#########6#########7#########8

rule all:
    input:
        #expand(config["data_dir"] + "/fastq/processed/{library_id}_proc_R1.fastq.gz", library_id = LIBRARY_IDS),
        #expand(config["data_dir"] + "/fastq/unpaired/{library_id}_unpr_R1.fastq.gz", library_id = LIBRARY_IDS),
        #expand(config["data_dir"] + "/fastq/processed/{library_id}_proc_R2.fastq.gz", library_id = LIBRARY_IDS),
        #expand(config["data_dir"] + "/fastq/unpaired/{library_id}_unpr_R2.fastq.gz", library_id = LIBRARY_IDS),
        #expand(config["data_dir"] + "/qc/{library_id}_{read}_fastqc.html", library_id = LIBRARY_IDS, read = ["R1","R2"]),
        #expand(config["data_dir"] + "/qc/{library_id}_proc_{read}_fastqc.html", library_id = LIBRARY_IDS, read = ["R1","R2"]),
        #expand(config["data_dir"] + "/bam/raw/{library_id}.bam", library_id = LIBRARY_IDS),
        #expand(config["data_dir"] + "/bam/raw/{library_id}.bam.bai", library_id = LIBRARY_IDS),
        #expand(config["data_dir"] + "/bam/filt/{library_id}_filt.bam",	library_id = LIBRARY_IDS),
        #expand(config["data_dir"] + "/bam/filt/{library_id}_filt.bam.bai", library_id = LIBRARY_IDS),
        #expand(config["data_dir"] + "/qc/{library_id}_samstats.txt", library_id = LIBRARY_IDS),
        #expand(config["data_dir"] + "/qc/{library_id}_flagstat.txt", library_id = LIBRARY_IDS),
        #expand(config["data_dir"] + "/qc/{library_id}_collect_wgs_metrics.txt", library_id = LIBRARY_IDS),
        #expand(config["data_dir"] + "/qc/{library_id}_deeptools_frag_lengths.txt", library_id = LIBRARY_IDS),
        #config["data_dir"] + "/qc/all_frag.tsv",
        #
        config["data_dir"] + "/qc/read_qc.tsv",
        get_ds_candidates,

rule symlink:
    input:
        config["data_dir"] + "/inputs/{library_id}_{read}.fastq.gz",
    output:
        config["data_dir"] + "/fastq/raw/{library_id}_{read}.fastq.gz",
    log:
        config["data_dir"] + "/logs/{library_id}_{read}_symlink.log"
    shell:
        """
        ln --force --relative --symbolic {input} {output} 2>{log}
        """

include: "read_preprocess.smk"
