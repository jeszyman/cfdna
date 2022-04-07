container: config["container"]

FQ_ID_RAW, = glob_wildcards(config["raw_fq_dir"] + "/{id}_R1.fastq.gz")

rule all:
    input:
        expand(config["processed_fq_dir"] + "/{fq_id}_proc_R1.fastq.gz", fq_id=FQ_ID_RAW),
        expand(config["processed_fq_dir"] + "/{fq_id}_proc_R2.fastq.gz", fq_id=FQ_ID_RAW),
        expand(config["unpr_fq_dir"] + "/{fq_id}_unpr_R1.fastq.gz", fq_id=FQ_ID_RAW),
        expand(config["unpr_fq_dir"] + "/{fq_id}_unpr_R2.fastq.gz", fq_id=FQ_ID_RAW),
        expand(config["qc_dir"] + "/{fq_id}_{read}_fastqc.html", fq_id=FQ_ID_RAW, read=["R1", "R2"]),
        expand(config["qc_dir"] + "/{fq_id}_proc_{read}_fastqc.html", fq_id=FQ_ID_RAW, read=["R1", "R2"]),
        expand(config["bam_dir"] + "/{fq_id}.sam", fq_id=FQ_ID_RAW),
        expand(config["bam_dir"] + "/{fq_id}.bam", fq_id=FQ_ID_RAW),
        expand(config["bam_dir"] + "/{fq_id}_dedup.bam", fq_id=FQ_ID_RAW),
        expand(config["bam_dir"] + "/{fq_id}_dedup_sort.bam", fq_id=FQ_ID_RAW),
        expand(config["bam_dir"] + "/{fq_id}_dedup_sort.bam.bai", fq_id=FQ_ID_RAW),

rule trimmomatic:
    input:
        read1 = config["raw_fq_dir"] + "/{fq_id}_R1.fastq.gz",
        read2 = config["raw_fq_dir"] + "/{fq_id}_R2.fastq.gz",
    params:
        adapter_fasta = config["inputs_dir"] + "/TruSeq3-PE.fa",
    output:
        read1 = config["processed_fq_dir"] + "/{fq_id}_proc_R1.fastq.gz",
        read1_unpr = config["unpr_fq_dir"] + "/{fq_id}_unpr_R1.fastq.gz",
        read2 = config["processed_fq_dir"] + "/{fq_id}_proc_R2.fastq.gz",
        read2_unpr = config["unpr_fq_dir"] + "/{fq_id}_unpr_R2.fastq.gz",	
    log:
        int = config["log_dir"] + "/trimmomatic_trimlog_{fq_id}.log",
        main = config["log_dir"] + "/trimmomatic_{fq_id}.log",
    shell:
        """
        workflow/scripts/trimmomatic.sh \
        {config[threads]} \
        {log.int} \
        {input.read1} \
        {input.read2} \
        {output.read1} \
        {output.read1_unpr} \
        {output.read2} \
        {output.read2_unpr} \
        {params.adapter_fasta} &> {log.main}
        """

rule fastqc_raw:
    input: 
        config["raw_fq_dir"] + "/{fq_id}_{read}.fastq.gz",
    params: 
        out_dir = config["qc_dir"],
    output:
        html = config["qc_dir"] + "/{fq_id}_{read}_fastqc.html", 
        zip = config["qc_dir"] + "/{fq_id}_{read}_fastqc.zip", 
    log: 
        config["log_dir"] + "/fastqc_raw_{fq_id}_{read}.log",
    shell:
        """
	fastqc --outdir {params.out_dir} \
	--quiet \
	--threads {config[threads]} {input} &> {log}
        """

rule fastqc_proc:
    input: 
        config["processed_fq_dir"] + "/{fq_id}_proc_{read}.fastq.gz",
    params: 
        out_dir = config["qc_dir"],
    output:
        html = config["qc_dir"] + "/{fq_id}_proc_{read}_fastqc.html", 
        zip = config["qc_dir"] + "/{fq_id}_proc_{read}_fastqc.zip", 
    log: 
        config["log_dir"] + "/fastqc_raw_{fq_id}_{read}.log",
    shell:
        """
	fastqc --outdir {params.out_dir} \
	--quiet \
	--threads {config[threads]} {input} &> {log}
        """

rule align:
    input:
        read1 = config["processed_fq_dir"] + "/{fq_id}_proc_R1.fastq.gz",
        read2 = config["processed_fq_dir"] + "/{fq_id}_proc_R2.fastq.gz",
    output:
        config["bam_dir"] + "/{fq_id}.sam",
    log:
        config["log_dir"] + "/align_{fq_id}.log"
    shell:
        """
        bwa mem -M -t 4 {config[bwa_index]} {input.read1} {input.read2} > {output}
	"""

rule alignment_processing:
    input:
        config["bam_dir"] + "/{fq_id}.sam",
    output:
        bam = config["bam_dir"] + "/{fq_id}.bam",
        dedup = config["bam_dir"] + "/{fq_id}_dedup.bam",
        sort = config["bam_dir"] + "/{fq_id}_dedup_sort.bam",
        index = config["bam_dir"] + "/{fq_id}_dedup_sort.bam.bai",
    log:
        config["log_dir"] + "/alignment_processing_{fq_id}.log"
    shell:
        """
        sambamba view -t {config[threads]} -S -f bam {input} > {output.bam}
        sambamba markdup -r -t {config[threads]} {output.bam} {output.dedup}
        sambamba sort -t {config[threads]} {output.dedup} -o {output.sort}
        sambamba index -t {config[threads]} {output.sort}
        """
