##################################################################
###   Integration testing snakefile for WGS cfDNA Processing   ###
##################################################################

import pandas as pd
import re
import numpy as np
container: config["container"]

# Setup sample name index as a python dictionary

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

LIBRARIES = list(lib_dict.keys())
FASTQS = list(lib_dict.values())

# List of downsampling values in millions of reads
MILREADS = config["MILREADS"]

# Makes the name bwa index directory from the config genome fasta
#  e.g. test/inputs/chr8.fa will make test/ref/chr8
genome_ref = config["genome_fasta"]
genome_ref = re.sub("inputs", lambda x: 'ref', genome_ref)
genome_ref = re.sub("\..*$", lambda x: '', genome_ref)

# Directory structure under data_dir:
cfdna_wgs_fastq_dir = config["data_dir"] + "/fastq/cfdna_wgs"
cfdna_wgs_bam_dir = config["data_dir"] + "/bam/cfdna_wgs"
cfdna_wgs_qc_dir = config["data_dir"] + "/qc/cfdna_wgs"
cfdna_wgs_log_dir = config["data_dir"] + "/logs/cfdna_wgs"

# Function acts on read_qc, generated in the workflow, to select libraries for
# downsampling. Notice library 2 does not downsample because it already has
# fewer than 3000 reads. Best practice for real data would be to use the
# MILREADS value in lieu of a specified number here.

def get_ds_candidates(wildcards):
    read_qc = pd.read_table(checkpoints.make_qc_tbl.get().output[0])
    test=read_qc.library[read_qc.reads_properly_paired > 3000].tolist()
    return expand(
	cfdna_wgs_bam_dir + "/ds/{library_id}_ds{milreads}.bam",
        library_id=test, milreads = MILREADS)

#########1#########2#########3#########4#########5#########6#########7#########8

rule all:
    input:
        #expand(cfdna_wgs_fastq_dir + "/raw/{library}_{read}.fastq.gz", library = lib_dict.keys(), read = ["R1", "R2"]),
        #expand(cfdna_wgs_fastq_dir + "/processed/{library_id}_proc_R1.fastq.gz", library_id = LIBRARIES),
        #expand(cfdna_wgs_fastq_dir + "/unpaired/{library_id}_unpr_R1.fastq.gz", library_id = LIBRARIES),
        #expand(cfdna_wgs_fastq_dir + "/processed/{library_id}_proc_R2.fastq.gz", library_id = LIBRARIES),
        #expand(cfdna_wgs_fastq_dir + "/unpaired/{library_id}_unpr_R2.fastq.gz", library_id = LIBRARIES),
        #expand(cfdna_wgs_qc_dir + "/{library_id}_{read}_fastqc.html", library_id = LIBRARIES, read = ["R1","R2"]),
        #expand(cfdna_wgs_qc_dir + "/{library_id}_proc_{read}_fastqc.html", library_id = LIBRARIES, read = ["R1","R2"]),
        #expand(cfdna_wgs_bam_dir + "/raw/{library_id}.bam", library_id = LIBRARIES),
        #expand(cfdna_wgs_bam_dir + "/raw/{library_id}.bam.bai", library_id = LIBRARIES),
        #expand(cfdna_wgs_bam_dir + "/filt/{library_id}_filt.bam", library_id = LIBRARIES),
        #expand(cfdna_wgs_bam_dir + "/filt/{library_id}_filt.bam.bai", library_id = LIBRARIES),
        #expand(cfdna_wgs_qc_dir + "/{library_id}_samstats.txt", library_id = LIBRARIES),
        #expand(cfdna_wgs_qc_dir + "/{library_id}_flagstat.txt", library_id = LIBRARIES),
        #expand(cfdna_wgs_qc_dir + "/{library_id}_collect_wgs_metrics.txt", library_id = LIBRARIES),
        #expand(cfdna_wgs_qc_dir + "/{library_id}_deeptools_frag_lengths.txt", library_id = LIBRARIES),
        #cfdna_wgs_qc_dir + "/all_frag.tsv",
        #
        # Final rules:
        cfdna_wgs_qc_dir + "/read_qc.tsv",
        get_ds_candidates,

rule symlink_inputs:
    input:
        lambda wildcards: lib_dict[wildcards.library],
    output:
        r1 = cfdna_wgs_fastq_dir + "/raw/{library}_R1.fastq.gz",
        r2 = cfdna_wgs_fastq_dir + "/raw/{library}_R2.fastq.gz",
    container:
        config["cfdna_wgs_container"]
    shell:
        """
        r2=$(echo {input} | sed "s/_R1/_R2/g")
        ln -sf --relative {input} {output.r1}
        ln -sf --relative $r2 {output.r2}
        """

include: "read_preprocess.smk"
