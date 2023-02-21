# [[file:~/repos/cfdna-wgs/cfdna-wgs.org::*Preamble][Preamble:1]]
##################################################################
###   Integration testing snakefile for WGS cfDNA Processing   ###
##################################################################
# Preamble:1 ends here

# [[file:~/repos/cfdna-wgs/cfdna-wgs.org::*Load%20packages][Load packages:1]]
import pandas as pd
import re
import numpy as np
# Load packages:1 ends here

# [[file:~/repos/cfdna-wgs/cfdna-wgs.org::*Variable%20naming][Variable naming:1]]
# Values directly from configuration file
threads = config["threads"]
FRAG_DISTROS = config["frag_distro"]
frag_threads = config["threads"]
genome_fasta = config["genome_fasta"]
frag_repo = config["frag_repo"]

# Directory values derived from data_dir in configuration YAML
data_dir                   = config["data_dir"]
frag                 = data_dir + "/analysis/frag"
frag_bams            = data_dir + "/analysis/frag/bams"
frag_fastqs          = data_dir + "/analysis/frag/fastqs"
frag_frag            = data_dir + "/analysis/frag/frag"
frag_frag_beds       = data_dir + "/analysis/frag/frag/beds"
frag_frag_counts     = data_dir + "/analysis/frag/frag/counts"
frag_frag_gc_distros = data_dir + "/analysis/frag/frag/distros"
qcdir                     = data_dir + "/analysis/qc"
benchdir                  = data_dir + "/benchmark"
logdir                    = data_dir + "/logs"
refdir                    = data_dir + "/re"

frag_scriptdir = config["frag_repo"] +  "/scripts"

bwa_dir = "{data_dir}/ref/hg38"
fasta_base = "GCA_000001405.15_GRCh38_no_alt_analysis_set"
frag_script_dir = "{frag_repo}/scripts"
# Variable naming:1 ends here

# [[file:~/repos/cfdna-wgs/cfdna-wgs.org::*Functions,%20miscellaneous][Functions, miscellaneous:1]]
#####################
###   Functions   ###
#####################

# Setup fragment sample name index as a python dictionary
frag_libs = pd.read_table("{data_dir}/test/inputs/frag_libs.tsv")

lib_path = "{data_dir}/test/inputs"

# Ensure readable fastqs
readable = []
for x in lib_path + "/" + frag_libs["r1_basename"]:
    readable.append(os.access(x, os.R_OK))
frag_libs['readable']=readable
frag_libs = frag_libs[frag_libs.readable == True]

# Make the dictionary
FRAG_LIBS = frag_libs["library"].tolist()
frag_libs_file_indict = lib_path + "/" + frag_libs["r1_basename"]
frag_lib_dict = dict(zip(FRAG_LIBS, frag_libs_file_indict))

# Make  a list of healthy libraries
FRAG_HEALTHY_LIBRARIES = frag_libs[frag_libs['cohort'] == 'healthy']['library'].tolist()
# Functions, miscellaneous:1 ends here

# [[file:~/repos/cfdna-wgs/cfdna-wgs.org::*All%20rule][All rule:1]]
rule all:
    input:
        expand("{data_dir}/analysis/frag/fastqs/{{library}}_raw_{{read}}.fastq.gz",
               library = list(frag_lib_dict.keys()),
               read = ["R1", "R2"]),
        "{data_dir}/ref/{fasta_base}.sa",
        #"{data_dir}/ref/{fasta_base}.sa",
        #logdir + "/aggregate_output",
        #frag_frag + "/ratios.tsv",
        #qcdir + "/frag_read_qc.tsv",
        #qcdir + "/frag_frag_len.tsv",
# All rule:1 ends here

# [[file:~/repos/cfdna-wgs/cfdna-wgs.org::*Benchmark%20aggregation][Benchmark aggregation:1]]
onsuccess:
    shell("""
        bash {frag_scriptdir}/agg_bench.sh {benchdir} {qcdir}/agg_bench.tsv
        """)
# Benchmark aggregation:1 ends here

# [[file:~/repos/cfdna-wgs/cfdna-wgs.org::*Symlink%20input%20fastqs][Symlink input fastqs:1]]
rule symlink_inputs:
    input:
        lambda wildcards: frag_lib_dict[wildcards.library],
    output:
        read1 = "{data_dir}/analysis/frag/fastqs/{{library}}_raw_R1.fastq.gz",
        read2 = "{data_dir}/analysis/frag/fastqs/{{library}}_raw_R2.fastq.gz",
    params:
        outdir = frag_fastqs,
        script = "{frag_script_dir}/symlink.sh",
    shell:
        """
        {params.script} \
        {input} \
        {output.read1} \
        {output.read2} \
        {params.outdir}
        """
# Symlink input fastqs:1 ends here

# [[file:~/repos/cfdna-wgs/cfdna-wgs.org::*Includes%20statements][Includes statements:1]]
include: frag_repo + "/workflow/frag_reads.smk"
#include: frag_repo + "/workflow/frag.smk"
# Includes statements:1 ends here
