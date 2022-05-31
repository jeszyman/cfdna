import pandas as pd

container: config["container"]

libraries = pd.read_table(config["data_dir"] + "/inputs/libraries.tsv")

LIBRARY_IDS = list(libraries.library.unique())

MILREADS = config["MILREADS"]

rule all:
    input:
        expand(config["data_dir"] + "/fastq/raw/{library_id}_{read}.fastq.gz", library_id = LIBRARY_IDS, read = ["R1", "R2"]),
        expand(config["data_dir"] + "/fastq/processed/{library_id}_proc_{read}.fastq.gz", library_id = LIBRARY_IDS, read = ["R1","R2"]),
        expand(config["data_dir"] + "/fastq/unpaired/{library_id}_unpr_R1.fastq.gz", library_id = LIBRARY_IDS, read = ["R1","R2"]),
        expand(config["data_dir"] + "/bam/{library_id}.sam", library_id = LIBRARY_IDS),
        expand(config["data_dir"] + "/qc/{library_id}_{read}_fastqc.html", library_id = LIBRARY_IDS, read = ["R1","R2"]),
        expand(config["data_dir"] + "/qc/{library_id}_proc_{read}_fastqc.html", library_id = LIBRARY_IDS, read = ["R1","R2"]),
        expand(config["data_dir"] + "/bam/{library_id}_dedup.bam", library_id = LIBRARY_IDS),
        expand(config["data_dir"] + "/bam/{library_id}_dedup.bam.bai", library_id = LIBRARY_IDS),
        expand(config["data_dir"] + "/qc/{library_id}_{bam_step}_samstats.txt", library_id = LIBRARY_IDS, bam_step= ["dedup","raw"]),
        expand(config["data_dir"] + "/qc/{library_id}_{bam_step}_flagstat.txt", library_id = LIBRARY_IDS, bam_step =["dedup","raw"]),
        expand(config["data_dir"] + "/bam/{library_id}_ds{milreads}.bam", library_id = LIBRARY_IDS, milreads = MILREADS),
        config["data_dir"] + "/qc/all_qc.html",

rule symlink:
    input:
        config["data_dir"] + "/inputs/{library_id}_{read}.fastq.gz",
    output:
        config["data_dir"] + "/fastq/raw/{library_id}_{read}.fastq.gz",
    shell:
        """
        ln --force --relative --symbolic {input} {output}
        """

include: "read_preprocess.smk"

rule multiqc:
    input:
        expand(config["data_dir"] + "/qc/{library_id}_{read}_fastqc.html", library_id = LIBRARY_IDS, read = ["R1","R2"]),
        expand(config["data_dir"] + "/qc/{library_id}_proc_{read}_fastqc.html", library_id = LIBRARY_IDS, read = ["R1","R2"]),
        expand(config["data_dir"] + "/qc/{library_id}_{bam_step}_samstats.txt", library_id = LIBRARY_IDS, bam_step= ["dedup","raw"]),
        expand(config["data_dir"] + "/qc/{library_id}_{bam_step}_flagstat.txt", library_id = LIBRARY_IDS, bam_step =["dedup","raw"]),
    params:
        out_dir = config["data_dir"] + "/qc"
    output:
        config["data_dir"] + "/qc/all_qc.html"
    shell:
        """
        multiqc {params.out_dir} \
        --force \
        --outdir {params.out_dir} \
        --filename all_qc
        """
