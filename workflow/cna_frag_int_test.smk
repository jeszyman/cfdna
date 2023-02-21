# [[file:~/repos/cfdna-wgs/cfdna-wgs.org::*Preamble][Preamble:1]]
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
# Preamble:1 ends here

# [[file:~/repos/cfdna-wgs/cfdna-wgs.org::*Variable%20naming][Variable naming:1]]
# Variable naming
benchdir = config["benchdir"]
frag_repo = config["frag_repo"]
frag_scriptdir = config["frag_scriptdir"]
logdir = config["logdir"]
threads = config["threads"]

# Suggested directory structure:
analysis = config["data_dir"] + "/analysis"
frag = config["data_dir"]      + "/analysis/frag"
frag_cna = config["data_dir"]  + "/analysis/frag/cna"
frag_frag = config["data_dir"] + "/analysis/frag/frag"

# Terminal variable paths:
#  (These variables are used directly in the cna snakefile)
frag_cna_in_bams      = frag_cna + "/input_bams"
frag_cna_frag_bams    = frag_cna + "/frag_bams"
frag_cna_wigs         = frag_cna + "/wigs"
frag_cna_ichor_nopon  = frag_cna + "/ichor_nopon"

frag_frag_input_bams  = frag_cna + "/input_bams"
frag_frag_beds       = frag_frag + "/beds"

frag_frag_counts     = frag_frag + "/counts"

refdir                 = config["data_dir"] + "/re"

# Additional variable names used directly in the cna snakefile:
chrom_sizes = config["chrom_sizes"]
genome_fasta = "/mnt/ris/aadel/Active/mpnst/inputs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"


#TMP_FRAG_LIBS = ["lib001_filt","lib002_filt"]

#chrs = "chr8"

chrs = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM",

keep_bed = refdir + "/hg38_keep.bed",
blklist = config["blklist"]
genome_ref = config["genome_re"]


FRAG_DISTROS = config["frag_distro"]

frag_threads = config["threads"]
frag_scriptdir = config["frag_scriptdir"]


frag_container = config["frag_container"]
default_container = config["default_container"]

autosome_bed = refdir + "/hg38_autosomes.bed",
frag_fastqs = frag + "/fastqs"
frag_bams = frag + "/bams"
qc = config["data_dir"] + "/qc"

# frag_container = config["frag_container"]


# frag_cna_bam_inputs   = config["dir"]["data"] + "/bam/filt"
# frag_cna_bam_fragfilt = config["dir"]["data"] + "/bam/frag"

# wig = config["dir"]["data"] + "/wig"
# ichor = config["dir"]["data"] + "/ichor"
# frag_logs = config["dir"]["data"] + "logs/frag"
# ichor_nopon = config["dir"]["data"] + "/ichor_nopon"
# Variable naming:1 ends here

# [[file:~/repos/cfdna-wgs/cfdna-wgs.org::*Functions][Functions:1]]
libraries = pd.read_table(config["data_dir"] + "/inputs/libraries.tsv")

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

FRAG_LIBS = list(lib_dict.keys())

cna_libraries = pd.read_table(config["data_dir"] + "/inputs/cna_libraries.tsv")

readable = []
for x in cna_libraries.bam_file:
    readable.append(os.access(x, os.R_OK))
cna_libraries['readable']=readable

cna_libraries = cna_libraries[cna_libraries.readable == True]

library_indict = cna_libraries["library"].tolist()
file_indict = cna_libraries["bam_file"].tolist()
lib_dict = dict(zip(library_indict, file_indict))

CNA_WGS_LIBRARIES = list(lib_dict.keys())
# Functions:1 ends here

# [[file:~/repos/cfdna-wgs/cfdna-wgs.org::*All%20rule][All rule:1]]
rule all:
    input:
# # From this snakefile:
#         # frag_symlink:
#         expand(frag_cna_in_bams +
#                "/{library}.bam",
#                library = lib_dict.keys()),
# # From cna.smk
#         # cna_frag_filt:
#         expand(frag_cna_frag_bams +
#                "/{library}_frag{frag_distro}.bam",
#                library = CNA_WGS_LIBRARIES,
#                frag_distro = FRAG_DISTROS),
#         # bam_to_wig:
#         expand(frag_cna_wigs +
#                "/{library}_frag{frag_distro}.wig",
#                library = CNA_WGS_LIBRARIES,
#                frag_distro = FRAG_DISTROS),
#         # ichor_nopon:
#         expand(frag_cna_ichor_nopon +
#                "/{library}_frag{frag_distro}.cna.seg",
#                library = CNA_WGS_LIBRARIES,
#                frag_distro = FRAG_DISTROS),
# From frag.smk
        # make_gc_map_bind:
        refdir + "/keep_5mb.bed",
        # filt_bam_to_frag_bed:
        expand(frag_frag_beds +
               "/{library}_filt.bed",
               library = CNA_WGS_LIBRARIES),
        # # gc_distro:
        # expand(frag_frag_gc_distros +
        #        "/{library}_gc_distro.csv",
        #        library = CNA_WGS_LIBRARIES),
        # # healthy_gc:
        # frag_frag_gc_distros + "/healthy_med.rds",
        # #
        # expand(frag_frag_beds +
        #        "/{library}_sampled_frag.bed",
        #       library = CNA_WGS_LIBRARIES),
        # expand(frag_frag_beds) /
        #        "{library}_norm_{length}.bed",
        #        library = CNA_WGS_LIBRARIES,
        #        length = ["short", "long"]),
        expand(frag_frag_counts +
               "/{library}_cnt_{length}.tmp",
               library = CNA_WGS_LIBRARIES,
               length = ["short", "long"]),
        frag_frag + "/frag_counts.tsv",
        #
        # unit_cent_sd:
        frag_frag + "/ratios.tsv",
# All rule:1 ends here

# [[file:~/repos/cfdna-wgs/cfdna-wgs.org::*Symlink%20input%20bams][Symlink input bams:1]]
# Symlink input bams
rule frag_symlink:
    container: frag_container,
    input: lambda wildcards: lib_dict[wildcards.library],
    output: frag_cna_in_bams + "/{library}.bam",
    shell:
        """
        ln --force --relative --symbolic {input} {output}
        """
# Symlink input bams:1 ends here

# [[file:~/repos/cfdna-wgs/cfdna-wgs.org::*Includes%20statements][Includes statements:1]]
include: frag_repo + "/workflow/reads.smk"
include: frag_repo + "/workflow/cna.smk"
include: frag_repo + "/workflow/frag.smk"
# Includes statements:1 ends here

# [[file:~/repos/cfdna-wgs/cfdna-wgs.org::*Filter%20downsampled%20bams%20to%20set%20fragment%20length%20distributions][Filter downsampled bams to set fragment length distributions:1]]
rule frag_filt:
    input:
        main = frag_bams + "/{library}_ds{downsample}.bam",
        check = logdir + "/{library}_{downsample}_made",
    output:
        nohead = temp(frag_bams + "/{library}_ds{downsample}_frag{frag_distro}.nohead"),
        onlyhead = temp(frag_bams + "/{library}_ds{downsample}_frag{frag_distro}.only"),
        final = frag_bams + "/{library}_ds{downsample}_frag{frag_distro}.bam",
    params:
        script = "{frag_script_dir}/frag_filt.sh",
        threads = frag_threads,
    shell:
        """
        frag_min=$(echo {wildcards.frag_distro} | sed -e "s/_.*$//g")
        frag_max=$(echo {wildcards.frag_distro} | sed -e "s/^.*_//g")
        {params.script} \
        {input.main} \
        {output.nohead} \
        $frag_min \
        $frag_max \
        {config[threads]} \
        {output.onlyhead} \
        {output.final}
        """
# Filter downsampled bams to set fragment length distributions:1 ends here

# [[file:~/repos/cfdna-wgs/cfdna-wgs.org::*Setup%20conditional%20execution%20of%20downsampled%20bams][Setup conditional execution of downsampled bams:1]]
# If downsample occured, then write filename into this per-library log, else leave the log file blank
rule log_dowsample:
    input: logdir + "/{library}_{downsample}_downsample.done",
    output: logdir + "/{library}_{downsample}_made",
    params:
        bamdir = frag_bams,
    shell:
        """
        dspath={params.bamdir}/{wildcards.library}_ds{wildcards.downsample}.bam
        if [ -f $dspath ]; then echo "$dspath"  > {output}; else touch {output}; fi
        """

# Use the downsampled bam logs to make a single text file of conditionally executed final targets.
# Specifically in this example, log text lines are in the form
# frag_bams + "/{library}_ds{downsample}_frag90_150.bam" to setup conditional execution of fragment filtering ONLY on downsampled bams
# Note alternative delimiter "~" to sed allows frag_wigs as param

checkpoint ds_cond_target_list:
    input: expand(logdir + "/{library}_{downsample}_made", library = FRAG_LIBS, downsample = DOWNSAMPLE),
    output: logdir + "/ds_final_targets",
    params:
        outdir = frag_bams,
        frag_distro=config["frag_distro"]
    shell:
        """
        if [ -f {output} ]; then rm {output}; fi
        cat {input} > {output}
        sed -i 's~^.*lib~{params.outdir}/lib~g' {output}
        sed -i 's/.bam$/_frag{params.frag_distro}.bam/g' {output}
        """

# Function jsut pulls the final target names out of ds_final_targets
def get_ds_targets(wildcards):
    with open(checkpoints.ds_cond_target_list.get(**wildcards).output[0], "r") as f:
      non_empty_files = [l.strip() for l in f.readlines()]
    return non_empty_files

# This rule allows execution of rules which will generate the conditional targets in ds_cond_target_list
rule make_ds_targets:
    input:
        get_ds_targets
    output: logdir + "/aggregate_output"
    run:
        with open(output[0], "w") as f:
            f.write("\n".join(input))
# Setup conditional execution of downsampled bams:1 ends here

# [[file:~/repos/cfdna-wgs/cfdna-wgs.org::*Downsample%20bams][Downsample bams:1]]
rule downsample_bams:
    input: frag_bams + "/{library}_filt.bam",
    output: touch(logdir + "/{library}_{downsample}_downsample.done"),
    params:
        out_dir = frag_bams,
        script = "{frag_script_dir}/downsample_bams.sh",
        suffix = "_filt.bam",
        threads = frag_threads,
    shell:
        """
        {params.script} \
        {input} \
        {wildcards.downsample} \
        {params.out_dir} \
        {params.suffix} \
        {params.threads}
        """
# Downsample bams:1 ends here
