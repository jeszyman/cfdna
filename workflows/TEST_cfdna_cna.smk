import pandas as pd
import os
from tabulate import tabulate

configfile: "./config/TEST_cfdna_cna.yaml"


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

class SampleTable:

    def __init__(self, tsv_path, selected_ids):
        # Load and filter the table to only include selected merged library IDs
        df = pd.read_csv(tsv_path, sep="\t")
        df = df[df["library_id"].isin(selected_ids)].copy()
        self.df = df

    @property
    def library_ids(self):
        return sorted(self.df["library_id"].unique())

    @property
    def frag_filt_ids(self):
        return sorted(f"{lib}.frag90_150" for lib in self.df["library_id"].unique())

#library_ids = ["kh_01"]

samples = SampleTable(config['input-tsv'], config['input-samples'])
print("Samples loaded:", samples.df["library_id"].unique().tolist())

rule all:
    input:
        # Filter cfDNA fragments to a length range (Rule cfdna_frag_filt)
        expand(f"{config['cfdna-cna-dir']}/bams/{{wkflow_id}}.frag{{frag_distro}}.bam",
               wkflow_id = samples.library_ids,
               frag_distro = ["90_150"]),
        #
        # Downsample length-filtered fragments (Rule cfdna_cna_downsample_bam)
        expand(f"{config['cfdna-cna-dir']}/ds-bams/{{wkflow_id}}.ds{{mil_reads}}.{{ext}}",
               wkflow_id = samples.library_ids,
               mil_reads = ["1"],
               ext = ["bam","bam.bai"]),

        expand(f"{config['cfdna-cna-dir']}/ds-bams/{{wkflow_id}}.ds{{mil_reads}}.{{ext}}",
               wkflow_id = samples.frag_filt_ids,
               mil_reads = ["1"],
               ext = [".bam",".bam.bai"]),

        #
        # Make wiggle file of coverage (Rule make_wig)
        expand(f"{config['cfdna-cna-dir']}/wigs/{{wkflow_id}}.frag{{frag_distro}}.ds{{mil_reads}}.wig",
               wkflow_id = samples.library_ids,
               frag_distro = ["90_150"],
               mil_reads = ["1"]),

        expand(f"{config['cfdna-cna-dir']}/wigs/{{wkflow_id}}.ds{{mil_reads}}.wig",
               wkflow_id = samples.library_ids,
               mil_reads = ["1"]),

        #
        # Run ichorCNA (Rule cfdna_ichor)
        expand(f"{config['cfdna-cna-dir']}/ichor/{{wkflow_id}}.frag{{frag_distro}}.ds{{mil_reads}}/{{wkflow_id}}.frag{{frag_distro}}.ds{{mil_reads}}.cna.seg",
               wkflow_id = samples.library_ids,
               frag_distro = ["90_150"],
               mil_reads = ["1"]),


include: f"{config['cfdna-repo-dir']}/workflows/cfdna_cna.smk"
