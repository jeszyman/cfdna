import pandas as pd
import os
from tabulate import tabulate

configfile: "./config/TEST_cfdna_wgs.yaml"

from pathlib import Path

def resolve_config_paths(config_dict):
    for k, v in config_dict.items():
        if isinstance(v, str):
            # Expand ~ and env vars
            config_dict[k] = os.path.expandvars(os.path.expanduser(v))
        elif isinstance(v, dict):
            resolve_config_paths(v)  # Recurse into nested dicts
        elif isinstance(v, list):
            config_dict[k] = [os.path.expandvars(os.path.expanduser(i)) if isinstance(i, str) else i for i in v]

resolve_config_paths(config)

LIBS = ["kh_01"]

cfdna_wgs_ref_names = ['ncbi_hg38']

rule all:
    input:
        # Process FASTQs
        expand(f"{config['data-dir']}/cfdna-wgs/fastqs/{{library_id}}.processed_{{read}}.fastq.gz",
               library_id = LIBS,
               read = ["R1","R2"]),

        # FastQC
        expand(f"{config['data-dir']}/cfdna-wgs/qc/{{library_id}}.{{processing}}_{{read}}_fastqc.{{sfx}}",
               library_id = LIBS,
               processing = ["raw", "processed"],
               read = ["R1", "R2"],
               sfx = ["html","zip"]),

        # BWA index
        expand(f"{config['data-dir']}/ref/bwa/{{name}}/{{name}}.amb",
               name = cfdna_wgs_ref_names),

        # Align
        expand(f"{config['cfdna-wgs-dir']}/bams/{{library_id}}.{{ref_name}}.bwa.coorsort.bam",
               library_id = LIBS,
               ref_name = cfdna_wgs_ref_names),

        # Post-process BAMs
        expand(f"{config['data-dir']}/cfdna-wgs/bams/{{library_id}}.{{ref_name}}.bwa.dedup.coorsort.{{suffix}}",
               library_id = LIBS,
               ref_name = cfdna_wgs_ref_names,
               suffix = ["bam", "bam.bai"]),

        expand(f"{config['data-dir']}/cfdna-wgs/bams/{{library_id}}.{{ref_name}}.bwa.dedup.coorsort.filt.{{suffix}}",
               library_id = LIBS,
               ref_name = cfdna_wgs_ref_names,
               suffix = ["bam", "bam.bai"]),

        # Alignment QC
        expand(f"{config['data-dir']}/cfdna-wgs/qc/{{library_id}}.{{ref_name}}.{{align_method}}.{{processing}}_{{stat}}.txt",
               library_id = LIBS,
               ref_name = cfdna_wgs_ref_names,
               align_method = ["bwa"],
               processing = ["coorsort","dedup.coorsort.filt"],
               stat = ["flagstat","samstat"])



include: f"{config['cfdna-repo']}/workflows/cfdna_wgs.smk"
