##################################################################
###   Integration testing snakefile for WGS cfDNA Processing   ###
##################################################################

import pandas as pd
import re
import numpy as np

# Suggested directory structure:
analysis      = config["datadir"] + "/analysis"
cfdna_wgs     = analysis          + "/cfdna_wgs"
logdir = config["datadir"] + "/logs"
refdir = config["datadir"] + "/ref"

# Terminal variable paths:
#  (These variables are used directly in the reads snakefile)
qc               = analysis + "/qc"
cfdna_wgs_fastqs = cfdna_wgs + "/fastqs"
cfdna_wgs_bams   = cfdna_wgs + "/bams"
autosome_bed = refdir + "/hg38_autosomes.bed",
keep_bed = refdir + "/hg38_keep.bed",

blacklist = config["blacklist"]
refdir = config["datadir"] + "/ref"
chrom_sizes = config["chrom_sizes"]
#chrs = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY",
chrs = "chr8"
FRAG_DISTROS = config["frag_distro"]
MILREADS = config["milreads"]
cfdna_wgs_threads = config["threads"]
cfdna_wgs_scriptdir = config["cfdna_wgs_scriptdir"]
default_container = config["default_container"]
cfdna_wgs_container = config["cfdna_wgs_container"]

# Makes the name bwa index directory from the config genome fasta
#  e.g. test/inputs/chr8.fa will make test/ref/chr8
genome_fasta = config["genome_fasta"]
genome_ref = config["genome_ref"]
blacklist = config["blacklist"]
genome_ref = config["genome_ref"]
# Directory structure under datadir:


cfdna_wgs_repo = config["cfdna_wgs_repo"]

# Setup sample name index as a python dictionary

libraries = pd.read_table(config["datadir"] + "/inputs/libraries.tsv")

readable = []
for x in libraries.file:
    readable.append(os.access(x, os.R_OK))
libraries['readable']=readable

cfdna_libraries = libraries
cfdna_libraries = cfdna_libraries[cfdna_libraries.library_type == "wgs"]
cfdna_libraries = cfdna_libraries[cfdna_libraries.isolation_type == "cfdna"]
cfdna_libraries = cfdna_libraries[cfdna_libraries.readable == True]

library_indict = cfdna_libraries["library"].tolist()
file_indict = cfdna_libraries["file"].tolist()
lib_dict = dict(zip(library_indict, file_indict))

CFDNA_WGS_LIBRARIES = list(lib_dict.keys())
CFDNA_WGS_FASTQS = list(lib_dict.values())


# Function acts on read_qc, generated in the workflow, to select libraries for
# downsampling.

def get_ds_candidates(wildcards):
    expand_milreads=1000000*MILREADS
    read_qc = pd.read_table(checkpoints.cfdna_wgs_make_qc_tsv.get().output[0])
    test = read_qc.library[read_qc.reads_mapped_and_paired_filt > expand_milreads].tolist()
    return expand(
        cfdna_wgs_bams + "/{library}_ds{milreads}.bam",
        library=test, milreads = MILREADS)

# Note, we could proceed with CNA if the resulting bams had enough reads using something like:
#POST_QC_LIBS, MILREADS_USED = glob_wildcards("test/analysis/cfdna_wgs_bams/{library}_ds0.002..bam")
#POST_QC_LIBS,= glob_wildcards("test/analysis/cfdna_wgs_bams/{library}_ds0.002.bam")

rule all:
    input:
# Read and alignment processing
        genome_ref,
        autosome_bed,
        keep_bed,
        expand(cfdna_wgs_fastqs +
               "/{library}_{processing}_{read}.fastq.gz",
               library = lib_dict.keys(),
               processing = ["raw", "processed", "unpaired"],
               read = ["R1", "R2"]),
        expand(cfdna_wgs_fastqs +
               "/{library}_failed_fastp.fastq.gz",
               library = CFDNA_WGS_LIBRARIES),
        expand(cfdna_wgs_bams +
               "/{library}_{processing}.bam",
               library = CFDNA_WGS_LIBRARIES,
               processing = ["raw", "dedup", "filt"]),
# Pipeline has a checkpoint here.
# Read and alignment QC
        expand(qc +
               "/{library}_{processing}_{read}_fastqc.html",
               library = CFDNA_WGS_LIBRARIES,
               processing = ["raw", "processed", "unpaired"],
               read = ["R1", "R2"]),
        expand(qc +
               "/{library}_{processing}_samstats.txt",
               library = CFDNA_WGS_LIBRARIES, processing = ["raw","dedup","filt"]),
        expand(qc +
               "/{library}_{processing}_flagstat.txt",
               library = CFDNA_WGS_LIBRARIES, processing = ["raw","dedup","filt"]),
        expand(qc +
               "/{library}_picard_depth.txt",
               library = CFDNA_WGS_LIBRARIES),
        qc + "/deeptools_frag_lengths.txt",
        qc + "/deeptools_frag_lengths.png",
        expand(qc +
               "/{library}_bamcoverage.bg",
               library = CFDNA_WGS_LIBRARIES),
        qc + "/cfdna_wgs_coverage.tsv",
        qc + "/cfdna_wgs_multiqc.html",
        qc + "/cfdna_wgs_read_qc.tsv",
        qc + "/cfdna_wgs_frag_len.tsv",
# Downsample bams
        #get_ds_candidates,

rule symlink_inputs:
    container: default_container,
    input:
        lambda wildcards: lib_dict[wildcards.library],
    output:
        read1 = cfdna_wgs_fastqs + "/{library}_raw_R1.fastq.gz",
        read2 = cfdna_wgs_fastqs + "/{library}_raw_R2.fastq.gz",
    params:
        outdir = cfdna_wgs_fastqs,
        script = cfdna_wgs_scriptdir + "/symlink.sh",
    shell:
        """
        {params.script} \
        {input} \
        {output.read1} \
        {output.read2} \
        {params.outdir}
        """

include: cfdna_wgs_repo + "/workflow/reads.smk"
