configfile: "./config/common.yaml"
READ_ID, = glob_wildcards(config["rename_dir"] + "/{id}_R1.fastq.gz")

rule all:
    input:
        expand(config["extracted_dir"] + "/{read_id}_R1.fastq", read_id=READ_ID),
        expand(config["extracted_dir"] + "/{read_id}_R2.fastq", read_id=READ_ID),
        expand(config["headfix_dir"] + "/{read_id}_{read}.fastq", read_id=READ_ID, read=["R1", "R2"]),

rule extract_cappseq_barcodes:
    input:
        read1 = config["rename_dir"] + "/{read_id}_R1.fastq.gz",
        read2 = config["rename_dir"] + "/{read_id}_R2.fastq.gz",
    output:
        read1 = config["rename_dir"] + "/{read_id}_R1.fastq",
        read2 = config["rename_dir"] + "/{read_id}_R2.fastq",	
        read1mv = config["extracted_dir"] + "/{read_id}_R1.fastq",
        read2mv = config["extracted_dir"] + "/{read_id}_R2.fastq",
    shell:
        """
        perl {config[cap_extract_script]} {input.read1} {input.read2}
        mv {output.read1} {output.read1mv}
        mv {output.read2} {output.read2mv}
        """

rule fix_headers:
    input:
	config["extracted_dir"] + "/{read_id}_{read}.fastq",
    output:
        unzip = config["headfix_dir"] + "/{read_id}_{read}.fastq",
	zip = config["headfix_dir"] + "/{read_id}_{read}.fastq.gz",
    shell:
        """
        cat {input} | awk '{if(NR%4==1){print substr($0, 1, length($0)-21)}else{print $0}}' > ${output.unzip}
        pigz -c -p {config[threads]} {output.unzip} > {output.zip} 
        """
