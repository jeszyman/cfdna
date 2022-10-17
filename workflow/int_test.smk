##################################################################
###   Integration testing snakefile for WGS cfDNA Processing   ###
##################################################################

import pandas as pd
import re
import numpy as np
#container: config["container"]

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

# Makes the name bwa index directory from the config genome fasta
#  e.g. test/inputs/chr8.fa will make test/ref/chr8
genome_ref = config["genome_fasta"]
genome_ref = re.sub("inputs", lambda x: 'ref', genome_ref)
genome_ref = re.sub("\..*$", lambda x: '', genome_ref)

# Directory structure under datadir:
cfdna_wgs_fastqs = config["datadir"] + "/cfdna-wgs-fastqs"
cfdna_wgs_bams =  config["datadir"] + "/cfdna-wgs-bams"
cfdna_wgs_failed_fastp = config["datadir"] + "/fastq_failed"

#########1#########2#########3#########4#########5#########6#########7#########8

rule all:
    input:
        expand(cfdna_wgs_fastqs +
           "/{library}_{processing}_{read}.fastq.gz",
            library = lib_dict.keys(),
            processing = ["raw", "processed", "unpaired"],
            read = ["R1", "R2"]),
        expand(cfdna_wgs_fastqs + "/{library}_failed_fastp.fastq.gz",
            library = CFDNA_WGS_LIBRARIES),
        expand(config["qcdir"] + "/{library}_{processing}_{read}_fastqc.html",
            library = CFDNA_WGS_LIBRARIES,
            processing = ["raw", "processed", "unpaired"],
            read = ["R1", "R2"]),
        genome_ref,
        expand(cfdna_wgs_bams + "/{library}_{processing}.bam",
            library = CFDNA_WGS_LIBRARIES,
            processing = ["raw", "dedup", "filt"]),
        expand(config["qcdir"] + "/{library}_{processing}_flagstat.txt",
            library = CFDNA_WGS_LIBRARIES,
            processing = ["raw", "dedup", "filt"]),
        expand(config["qcdir"] + "/{library}_{processing}_samstats.txt",
            library = CFDNA_WGS_LIBRARIES,
            processing = ["raw", "dedup", "filt"]),
        expand(config["qcdir"] + "/{library}_picard_wgs.txt", library = CFDNA_WGS_LIBRARIES),
        config["qcdir"] + "/deeptools_frag_lengths.png",
        config["qcdir"] + "/deeptools_frag_lengths.txt",
        expand(config["qcdir"] + "/{library}_bamcoverage.bg", library = CFDNA_WGS_LIBRARIES),
        config["qcdir"] + "/coverage.tsv",
        config["qcdir"] + "/coverage.pdf",
        config["qcdir"] + "/all_cfdna_wgs.html",
        config["qcdir"] + "/cfdna_wgs_read_qc.tsv",
        config["qcdir"] + "/cfdna_wgs_read_qc.tsv",
        config["qcdir"] + "/cfdna_wgs_frag_len.tsv",

rule symlink_inputs:
    container:
        config["container"]["default"],
    input:
        lambda wildcards: lib_dict[wildcards.library],
    output:
        read1 = cfdna_wgs_fastqs + "/{library}_raw_R1.fastq.gz",
        read2 = cfdna_wgs_fastqs + "/{library}_raw_R2.fastq.gz",
    params:
        outdir = cfdna_wgs_fastqs,
        script = config["scriptdir"]["cfdna_wgs"] + "/symlink.sh",
    shell:
        """
        {params.script} \
        {input} \
        {output.read1} \
        {output.read2} \
        {params.outdir}
        """

include: config["repo"]["cfdna_wgs"] + "/workflow/cfdna_wgs.smk"
