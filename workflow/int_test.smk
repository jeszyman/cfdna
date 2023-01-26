##################################################################
###   Integration testing snakefile for WGS cfDNA Processing   ###
##################################################################

import pandas as pd
import re
import numpy as np

# Values directly from configuration file
DOWNSAMPLE = config["downsample"]
threads = config["threads"]
FRAG_DISTROS = config["frag_distro"]
cfdna_wgs_threads = config["threads"]
default_container = config["default_container"]
cfdna_wgs_container = config["cfdna_wgs_container"]
genome_fasta = config["genome_fasta"]
genome_ref = config["genome_ref"]
cfdna_wgs_repo = config["cfdna_wgs_repo"]

# Directory values derived from datadir in configuration YAML
datadir                   = config["datadir"]
cfdna_wgs                 = datadir + "/analysis/cfdna_wgs"
cfdna_wgs_bams            = datadir + "/analysis/cfdna_wgs/bams"
cfdna_wgs_fastqs          = datadir + "/analysis/cfdna_wgs/fastqs"
cfdna_wgs_frag            = datadir + "/analysis/cfdna_wgs/frag"
cfdna_wgs_frag_beds       = datadir + "/analysis/cfdna_wgs/frag/beds"
cfdna_wgs_frag_counts     = datadir + "/analysis/cfdna_wgs/frag/counts"
cfdna_wgs_frag_gc_distros = datadir + "/analysis/cfdna_wgs/frag/distros"
qcdir                     = datadir + "/analysis/qc"
benchdir                  = datadir + "/benchmark"
logdir                    = datadir + "/logs"
refdir                    = datadir + "/ref"

cfdna_wgs_scriptdir = config["cfdna_wgs_repo"] +  "/scripts"

###   Functions   ###

# Setup sample name index as a python dictionary
cfdna_wgs_libraries = pd.read_table(config["datadir"] + "/inputs/libraries.tsv")

readable = []
for x in cfdna_wgs_libraries.file:
    readable.append(os.access(x, os.R_OK))
# Ensure readable fastqs
cfdna_wgs_libraries['readable']=readable
cfdna__wgs_libraries = cfdna_wgs_libraries[cfdna_wgs_libraries.readable == True]
# Ensure correct library type per sample sheet
cfdna_wgs_libraries = cfdna_wgs_libraries[cfdna_wgs_libraries.library_type == "wgs"]
cfdna_wgs_libraries = cfdna_wgs_libraries[cfdna_wgs_libraries.isolation_type == "cfdna"]

# Make the dictionary
cfdna_wgs_library_indict = cfdna_wgs_libraries["library"].tolist()
cfdna_wgs_file_indict = cfdna_wgs_libraries["file"].tolist()
cfdna_wgs_lib_dict = dict(zip(cfdna_wgs_library_indict, cfdna_wgs_file_indict))

CFDNA_WGS_LIBRARIES = list(cfdna_wgs_lib_dict.keys())
CFDNA_WGS_FASTQS = list(cfdna_wgs_lib_dict.values())

# Make  a list of healthy libraries
CFDNA_WGS_HEALTHY_LIBRARIES = cfdna_wgs_libraries[cfdna_wgs_libraries['cohort'] == 'healthy']['library'].tolist()

rule all:
    input:
        logdir + "/aggregate_output",
        cfdna_wgs_frag + "/ratios.tsv",
        qcdir + "/cfdna_wgs_read_qc.tsv",
        qcdir + "/cfdna_wgs_frag_len.tsv",

onsuccess:
    shell("""
        bash {cfdna_wgs_scriptdir}/agg_bench.sh {benchdir} {qcdir}/agg_bench.tsv
        """)

rule symlink_inputs:
    container: default_container,
    input:
        lambda wildcards: cfdna_wgs_lib_dict[wildcards.library],
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
