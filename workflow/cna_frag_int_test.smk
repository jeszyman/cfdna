#########1#########2#########3#########4#########5#########6#########7#########8
#                                                                              #
#      Integration Testing Snakefile for Analysis of Cell-free DNA             #
#    Whole Genome Sequencing Copy Number Alteration and Fragmentomics          #
#                                                                              #
#########1#########2#########3#########4#########5#########6#########7#########8

# Load necessary packages for snakemake run
import pandas as pd
import re
import numpy as np

# Variable naming
benchdir = config["benchdir"]
cfdna_wgs_repo = config["cfdna_wgs_repo"]
cfdna_wgs_scriptdir = config["cfdna_wgs_scriptdir"]
logdir = config["logdir"]
threads = config["threads"]

# Suggested directory structure:
analysis = config["datadir"] + "/analysis"
cfdna_wgs = config["datadir"]      + "/analysis/cfdna_wgs"
cfdna_wgs_cna = config["datadir"]  + "/analysis/cfdna_wgs/cna"
cfdna_wgs_frag = config["datadir"] + "/analysis/cfdna_wgs/frag"

# Terminal variable paths:
#  (These variables are used directly in the cna snakefile)
cfdna_wgs_cna_in_bams      = cfdna_wgs_cna + "/input_bams"
cfdna_wgs_cna_frag_bams    = cfdna_wgs_cna + "/frag_bams"
cfdna_wgs_cna_wigs         = cfdna_wgs_cna + "/wigs"
cfdna_wgs_cna_ichor_nopon  = cfdna_wgs_cna + "/ichor_nopon"

cfdna_wgs_frag_input_bams  = cfdna_wgs_cna + "/input_bams"
cfdna_wgs_frag_beds       = cfdna_wgs_frag + "/beds"
cfdna_wgs_frag_gc_distros = cfdna_wgs_frag + "/distros"
cfdna_wgs_frag_counts     = cfdna_wgs_frag + "/counts"

refdir                 = config["datadir"] + "/ref"

# Additional variable names used directly in the cna snakefile:
chrom_sizes = config["chrom_sizes"]
genome_fasta = "/mnt/ris/aadel/Active/mpnst/inputs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

CFDNA_WGS_HEALTHY_LIBRARIES = ["lib003", "lib004"]
#TMP_FRAG_LIBS = ["lib001_filt","lib002_filt"]

#chrs = "chr8"

chrs = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM",

keep_bed = refdir + "/hg38_keep.bed",
blklist = config["blklist"]
genome_ref = config["genome_ref"]


FRAG_DISTROS = config["frag_distro"]

cfdna_wgs_threads = config["threads"]
cfdna_wgs_scriptdir = config["cfdna_wgs_scriptdir"]

default_container = config["default_container"]
cfdna_wgs_container = config["cfdna_wgs_container"]


autosome_bed = refdir + "/hg38_autosomes.bed",
cfdna_wgs_fastqs = cfdna_wgs + "/fastqs"
cfdna_wgs_bams = cfdna_wgs + "/bams"
qc = config["datadir"] + "/qc"

# cfdna_wgs_container = config["cfdna_wgs_container"]


# cfdna_wgs_cna_bam_inputs   = config["dir"]["data"] + "/bam/filt"
# cfdna_wgs_cna_bam_fragfilt = config["dir"]["data"] + "/bam/frag"

# wig = config["dir"]["data"] + "/wig"
# ichor = config["dir"]["data"] + "/ichor"
# cfdna_wgs_logs = config["dir"]["data"] + "logs/cfdna_wgs"
# ichor_nopon = config["dir"]["data"] + "/ichor_nopon"

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

cna_libraries = pd.read_table(config["datadir"] + "/inputs/cna_libraries.tsv")

readable = []
for x in cna_libraries.bam_file:
    readable.append(os.access(x, os.R_OK))
cna_libraries['readable']=readable

cna_libraries = cna_libraries[cna_libraries.readable == True]

library_indict = cna_libraries["library"].tolist()
file_indict = cna_libraries["bam_file"].tolist()
lib_dict = dict(zip(library_indict, file_indict))

CNA_WGS_LIBRARIES = list(lib_dict.keys())

rule all:
    input:
# # From this snakefile:
#         # cfdna_wgs_symlink:
#         expand(cfdna_wgs_cna_in_bams +
#                "/{library}.bam",
#                library = lib_dict.keys()),
# # From cna.smk
#         # cna_frag_filt:
#         expand(cfdna_wgs_cna_frag_bams +
#                "/{library}_frag{frag_distro}.bam",
#                library = CNA_WGS_LIBRARIES,
#                frag_distro = FRAG_DISTROS),
#         # bam_to_wig:
#         expand(cfdna_wgs_cna_wigs +
#                "/{library}_frag{frag_distro}.wig",
#                library = CNA_WGS_LIBRARIES,
#                frag_distro = FRAG_DISTROS),
#         # ichor_nopon:
#         expand(cfdna_wgs_cna_ichor_nopon +
#                "/{library}_frag{frag_distro}.cna.seg",
#                library = CNA_WGS_LIBRARIES,
#                frag_distro = FRAG_DISTROS),
# From frag.smk
        # make_gc_map_bind:
        refdir + "/keep_5mb.bed",
        # filt_bam_to_frag_bed:
        expand(cfdna_wgs_frag_beds +
               "/{library}_filt.bed",
               library = CNA_WGS_LIBRARIES),
        # # gc_distro:
        # expand(cfdna_wgs_frag_gc_distros +
        #        "/{library}_gc_distro.csv",
        #        library = CNA_WGS_LIBRARIES),
        # # healthy_gc:
        # cfdna_wgs_frag_gc_distros + "/healthy_med.rds",
        # #
        # expand(cfdna_wgs_frag_beds +
        #        "/{library}_sampled_frag.bed",
        #       library = CNA_WGS_LIBRARIES),
        # expand(cfdna_wgs_frag_beds) /
        #        "{library}_norm_{length}.bed",
        #        library = CNA_WGS_LIBRARIES,
        #        length = ["short", "long"]),
        expand(cfdna_wgs_frag_counts +
               "/{library}_cnt_{length}.tmp",
               library = CNA_WGS_LIBRARIES,
               length = ["short", "long"]),
        cfdna_wgs_frag + "/frag_counts.tsv",
        #
        # unit_cent_sd:
        cfdna_wgs_frag + "/ratios.tsv",

# Symlink input bams
rule cfdna_wgs_symlink:
    container: cfdna_wgs_container,
    input: lambda wildcards: lib_dict[wildcards.library],
    output: cfdna_wgs_cna_in_bams + "/{library}.bam",
    shell:
        """
        ln --force --relative --symbolic {input} {output}
        """

include: cfdna_wgs_repo + "/workflow/reads.smk"
include: cfdna_wgs_repo + "/workflow/cna.smk"
include: cfdna_wgs_repo + "/workflow/frag.smk"
