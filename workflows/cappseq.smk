configfile: "./config/common.yaml"
demultiplex_fastqs, = glob_wildcards(config["demultiplex_dir"] + "/{id}.fastq.gz")

rule all:
    input: expand(config["rename_dir"] + "{id}.fastq.gz", id=demultiplex_fastqs)

rule rename:
    input: config["demultiplex_dir"] + "/{id}.fastq.gz",
    output: config["rename_dir"] + "/{id}.fastq.gz",
    shell:
        """
        ln -s {input} {output}
        """
