configfile: "./config/common.yaml"
READ_ID, = glob_wildcards(config["rename_dir"] + "/{id}_R1.fastq.gz")

rule all:
    input:
        expand(config["extracted_dir"] + "/{read_id}_extract_{read}.fastq.gz", read_id=READ_ID, read=["R1","R2"]),

rule extract_cappseq_barcodes:
    input:
        read1 = config["rename_dir"] + "/{read_id}_R1.fastq.gz",
        read2 = config["rename_dir"] + "/{read_id}_R2.fastq.gz",
    output:
        read1 = temp(config["rename_dir"] + "/{read_id}_R1.fastq"),
        read2 = temp(config["rename_dir"] + "/{read_id}_R2.fastq"),
        read1mv = config["extracted_dir"] + "/{read_id}_R1.fastq",
        read2mv = config["extracted_dir"] + "/{read_id}_R2.fastq",
    resources:
        mem_mb=10000	
    shell:
        """
        perl {config[cap_extract_script]} {input.read1} {input.read2}
        cp {output.read1} {output.read1mv}
        cp {output.read2} {output.read2mv}
        """

rule fix_headers:
    input:
        config["extracted_dir"] + "/{read_id}_{read}.fastq",
    output:
        unzip = config["extracted_dir"] + "/{read_id}_extract_{read}.fastq",
        zip = config["extracted_dir"] + "/{read_id}_extract_{read}.fastq.gz",	
    shell:
        """
        workflow/scripts/test.sh {input} {output.unzip}
        pigz -c -p {config[threads]} {output.unzip} > {output.zip} 
        """
