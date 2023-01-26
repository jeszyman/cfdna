##################################################################
###   Integration testing snakefile for WGS cfDNA Processing   ###
##################################################################

import pandas as pd
import re
import numpy as np

cfdna_wgs_scriptdir = config["cfdna_wgs_repo"] +  "/scripts"
CFDNA_WGS_HEALTHY_LIBRARIES = ["lib003", "lib004"]

# Values directly from configuration YAML
benchdir = config["benchdir"]
DOWNSAMPLE = config["downsample"]
datadir = config["datadir"]
qcdir = config["qcdir"]

# Suggested directory structure:
analysis      = config["datadir"] + "/analysis"
cfdna_wgs     = analysis          + "/cfdna_wgs"
logdir = config["datadir"] + "/logs"
refdir = config["datadir"] + "/ref"

# Terminal variable paths:
#  (These variables are used directly in the reads snakefile)
qc               = analysis + "/qc"
qcdir = config["datadir"] + "/analysis/qc"
cfdna_wgs_fastqs = cfdna_wgs + "/fastqs"
cfdna_wgs_bams   = cfdna_wgs + "/bams"
autosome_bed = refdir + "/hg38_autosomes.bed",
keep_bed = refdir + "/hg38_keep.bed",

cfdna_wgs_frag = cfdna_wgs + "/frag"
cfdna_wgs_frag_beds = cfdna_wgs_frag + "/beds"
cfdna_wgs_frag_gc_distros = cfdna_wgs_frag + "/distros"
cfdna_wgs_frag_counts = cfdna_wgs_frag + "/counts"

threads = config["threads"]


#TMP
cfdna_wgs_cna_in_bams = cfdna_wgs_bams
cfdna_wgs_cna_frag_bams = cfdna_wgs_bams
cfdna_wgs_cna_wigs = cfdna_wgs_bams
cfdna_wgs_cna_ichor_nopon = cfdna_wgs_bams

cfdna_wgs_wigs = cfdna_wgs + "/wigs"

cfdna_wgs_ichor_nopon = cfdna_wgs + "/ichor_nopon"
refdir = config["datadir"] + "/ref"

chrs = "chr8"
FRAG_DISTROS = config["frag_distro"]
cfdna_wgs_threads = config["threads"]

default_container = config["default_container"]
cfdna_wgs_container = config["cfdna_wgs_container"]



# Makes the name bwa index directory from the config genome fasta
#  e.g. test/inputs/chr8.fa will make test/ref/chr8
genome_fasta = config["genome_fasta"]
genome_ref = config["genome_ref"]

genome_ref = config["genome_ref"]
# Directory structure under datadir:

genome_fasta_cna = "/mnt/ris/aadel/Active/mpnst/inputs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

#chrs = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM",


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

#+begin_src snakemake
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

rule all:
    input:
        logdir + "/aggregate_output",
        cfdna_wgs_frag + "/ratios.tsv",
        qcdir + "/cfdna_wgs_read_qc.tsv",
        qcdir + "/cfdna_wgs_frag_len.tsv",

onsuccess:
    shell("""
        bash {cfdna_wgs_scriptdir}/agg_bench.sh {benchdir} {qc}/agg_bench.tsv
        """)

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
include: cfdna_wgs_repo + "/workflow/frag.smk"
