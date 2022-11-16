##################################################################
###   Integration testing snakefile for WGS cfDNA Processing   ###
##################################################################

import pandas as pd
import re
import numpy as np

blacklist = config["blacklist"]
refdir = config["datadir"] + "/ref"
chrom_sizes = config["chrom_sizes"]
#chrs = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY",
chrs = "chr8"
FRAG_DISTROS = config["frag_distro"]
MILREADS = config["milreads"]
cfdna_wgs_threads = config["threads"]
cfdna_wgs_scriptdir = config["cfdna_wgs_scriptdir"]
analysis = config["datadir"] + "/analysis"
default_container = config["default_container"]
cfdna_wgs_container = config["cfdna_wgs_container"]
logdir = config["datadir"] + "/logs"
# Makes the name bwa index directory from the config genome fasta
#  e.g. test/inputs/chr8.fa will make test/ref/chr8
genome_fasta = config["genome_fasta"]
genome_ref = config["genome_fasta"]
genome_ref = re.sub("inputs", lambda x: 'ref', genome_ref)
genome_ref = re.sub("\..*$", lambda x: '', genome_ref)

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
        analysis + "/cfdna_wgs_bams/{library}_ds{milreads}.bam",
        library=test, milreads = MILREADS)

#POST_QC_LIBS, MILREADS_USED = glob_wildcards("test/analysis/cfdna_wgs_bams/{library}_ds{milreads}.bam")

rule all:
    input:
        expand(analysis + "/cfdna_wgs_fastqs/{library}_{processing}_{read}.fastq.gz",
            library = lib_dict.keys(),
            processing = ["raw", "processed", "unpaired"],
            read = ["R1", "R2"]),
        expand(analysis + "/cfdna_wgs_fastqs/{library}_failed_fastp.fastq.gz",
             library = CFDNA_WGS_LIBRARIES),
        expand(analysis + "/qc/{library}_{processing}_{read}_fastqc.html",
             library = CFDNA_WGS_LIBRARIES,
             processing = ["raw", "processed", "unpaired"],
             read = ["R1", "R2"]),
        genome_ref,
        expand(analysis + "/cfdna_wgs_bams/{library}_{processing}.bam",
             library = CFDNA_WGS_LIBRARIES,
             processing = ["raw", "dedup", "filt"]),
        expand(analysis + "/qc/{library}_{processing}_samstats.txt",
               library = CFDNA_WGS_LIBRARIES, processing = ["raw","dedup","filt"]),
        expand(analysis + "/qc/{library}_{processing}_flagstat.txt",
               library = CFDNA_WGS_LIBRARIES, processing = ["raw","dedup","filt"]),
        expand(analysis + "/qc/{library}_picard_depth.txt", library = CFDNA_WGS_LIBRARIES),
        analysis + "/qc/deeptools_frag_lengths.txt",
        analysis + "/qc/deeptools_frag_lengths.png",
        expand(analysis + "/qc/{library}_bamcoverage.bg", library = CFDNA_WGS_LIBRARIES),
        analysis + "/qc/cfdna_wgs_coverage.tsv",
        analysis + "/qc/cfdna_wgs_multiqc.html",
        analysis + "/qc/cfdna_wgs_read_qc.tsv",
        analysis + "/qc/cfdna_wgs_frag_len.tsv",
        get_ds_candidates,
        # expand(analysis + "/cfdna_wgs_frag/{library}_ds{milreads}_frag{frag_distro}.bam",
        #     library = POST_QC_LIBS, milreads = MILREADS_USED, frag_distro = FRAG_DISTROS),
        # expand(analysis + "/cfdna_wgs_frag/{library}_ds{milreads}_frag{frag_distro}.wig",
        #     library = POST_QC_LIBS, milreads = MILREADS_USED, frag_distro = FRAG_DISTROS),

rule symlink_inputs:
    container: default_container,
    input:
        lambda wildcards: lib_dict[wildcards.library],
    output:
        read1 = analysis + "/cfdna_wgs_fastqs/{library}_raw_R1.fastq.gz",
        read2 = analysis + "/cfdna_wgs_fastqs/{library}_raw_R2.fastq.gz",
    params:
        outdir = analysis + "/cfdna_wgs_fastqs",
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
include: cfdna_wgs_repo + "/workflow/cna.smk"
