rule trimmomatic:
    input:
        read1 = config["fq_dir"] + "/{read_id}_R1.fastq.gz",
        read2 = config["fq_dir"] + "/{read_id}_R2.fastq.gz",
    params:
        adapter_fasta = config["inputs_dir"] + "/TruSeq3-PE.fa",
    output:
        read1 = config["processed_fq_dir"] + "/{read_id}_proc_R1.fastq.gz",
        read1_unpr = config["unpr_fq_dir"] + "/{read_id}_unpr_R1.fastq.gz",
        read2 = config["processed_fq_dir"] + "/{read_id}_proc_R2.fastq.gz",
        read2_unpr = config["unpr_fq_dir"] + "/{read_id}_unpr_R2.fastq.gz",	
    log:
        int = config["log_dir"] + "/trimmomatic_trimlog_{read_id}.log",
        main = config["log_dir"] + "/trimmomatic_{read_id}.log",
    shell:
        """
        trimmomatic PE \
                    -threads {config[threads]} \
                    -trimlog {log.int} \
                    {input.read1} {input.read2} \
                    {output.read1} {output.read1_unpr} \
                    {output.read2} {output.read2_unpr} \
                    ILLUMINACLIP:{params.adapter_fasta}:2:30:10 \
                    LEADING:10 TRAILING:10 MAXINFO:50:0.97 MINLEN:20 &> {log.main}
        """

rule align:
    input:
        read1 = config["processed_fq_dir"] + "/{read_id}_proc_R1.fastq.gz",
        read2 = config["processed_fq_dir"] + "/{read_id}_proc_R2.fastq.gz",
    output:
        config["bam_dir"] + "/{read_id}.sam",
    log:
        config["log_dir"] + "/align_{read_id}.log"
    shell:
        """
        bwa mem -M -t 4 {config[bwa_index]} {input.read1} {input.read2} > {output}
	"""

rule fastqc:
    input: 
        raw=config["fq_dir"] + "/{read_id}_{read}.fastq.gz",
        proc=config["processed_fq_dir"] + "/{read_id}_proc_{read}.fastq.gz",
    params: 
        out_dir = config["qc_dir"],
    output:
        raw_html = config["qc_dir"] + "/{read_id}_{read}_fastqc.html",
        proc_html = config["qc_dir"] + "/{read_id}_proc_{read}_fastqc.html", 	
    log: 
        raw = config["log_dir"] + "/fastqc_raw_{read_id}_{read}.log",
        proc = config["log_dir"] + "/fastqc_proc_{read_id}_{read}.log",	
    shell:
        """
        fastqc --outdir {params.out_dir} \
        --quiet \
        --threads {config[threads]} {input.raw} &> {log}
        fastqc --outdir {params.out_dir} \
        --quiet \
        --threads {config[threads]} {input.proc} &> {log}
        """

rule alignment_processing:
    input:
        config["bam_dir"] + "/{read_id}.sam",
    output:
        bam = config["bam_dir"] + "/{read_id}_raw.bam",
        dedup = temp(config["bam_dir"] + "/{read_id}_dedup_unsort.bam"),
        sort = config["bam_dir"] + "/{read_id}_dedup.bam",
        index = config["bam_dir"] + "/{read_id}_dedup.bam.bai",
    log:
        config["log_dir"] + "/alignment_processing_{read_id}.log"
    shell:
        """
        sambamba view -t {config[threads]} -S -f bam {input} > {output.bam}
        sambamba markdup -r -t {config[threads]} {output.bam} {output.dedup}
        sambamba sort -t {config[threads]} {output.dedup} -o {output.sort}
        sambamba index -t {config[threads]} {output.sort}
        """

rule alignment_qc:
    input:
        config["bam_dir"] + "/{read_id}_{bam_step}.bam",
    output:
        samstat = config["qc_dir"] + "/{read_id}_{bam_step}_samstats.txt",
        flagstat = config["qc_dir"] + "/{read_id}_{bam_step}_flagstat.txt",        
    shell:
        """
        samtools stats {input} > {output.samstat}
        samtools flagstat {input} > {output.flagstat}
        """

rule downsample_bams:
    input:
        bam = config["bam_dir"] + "/{read_id}_dedup.bam",
    output:
        config["bam_dir"] + "/{read_id}_ds{milreads}.bam",
    shell:
        """
        reads=$(echo {wildcards.milreads}000000)
        workflow/scripts/downsample_bam.sh {input} $reads {output}
        """
