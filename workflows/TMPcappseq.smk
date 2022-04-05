configfile: "./config/common.yaml"
demultiplex_fastqs, = glob_wildcards(config["demultiplex_dir"] + "/{id}.fastq.gz")
demultiplex_fastqs, = glob_wildcards(config["demultiplex_dir"] + "/{id_pt1}_{read}_{id_pt2}.fastq.gz")

rule all:
    input: expand(config["rename_dir"] + "/{id}.fastq.gz", id=demultiplex_fastqs)

rule copy:
    input: config["demultiplex_dir"] + "/{id}.fastq.gz",
    output: config["rename_dir"] + "/{id}.fastq.gz",
    shell:
        """
        cp {input} {output}
        """

rule rename1:
    input:
        read1 = config["rename_dir"] + "/{id}_R1_{bcode}.fastq.gz",
        read2 = config["rename_dir"] + "/{id}_R2_{bcode}.fastq.gz",	
    shell:
        """
        rename s/\.fastq.gz/_R1.fastq.gz/g {input.read1}
        rename s/\.fastq.gz/_R2.fastq.gz/g {input.read2}
        """

rule rename2:
    input:
        read1 = config["rename_dir"] + "/{id}_R1_{bcode}_R1.fastq.gz",
        read2 = config["rename_dir"] + "/{id}_R2_{bcode}_R2.fastq.gz",
    shell:
        """
        rename s/_R1_/_/g {input.read1}
        rename s/_R2_/_/g {input.read2}
        """

rule bar_extract:
    input:
        read1 = config

	
