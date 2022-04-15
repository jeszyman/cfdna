container: config["container"]
import pandas as pd
import numpy as np

samples = pd.read_table(config["inputs_dir"] + "/samples.tsv")
sampledict = dict(zip(samples['new_name'], samples['old_name']))

wildcard_constraints:
    read_id='|'.join([re.escape(x) for x in sampledict.keys()]),

rule all:
    input:
        expand(config["processed_fq_dir"] + "/{read_id}_proc_{read}.fastq.gz", read_id = sampledict.keys(), read = ["R1","R2"]),
        expand(config["unpr_fq_dir"] + "/{read_id}_unpr_R1.fastq.gz", read_id = sampledict.keys(), read = ["R1","R2"]),
	expand(config["bam_dir"] + "/{read_id}_dedup.bam", read_id = sampledict.keys()),
        config["qc_dir"] + "/all_qc.html",
        expand(config["bam_dir"] + "/{read_id}_ds{milreads}.bam", read_id = sampledict.keys(), milreads = config["MILREADS"]),

rule rename:
    params:
        old_sample_id=lambda wcs: sampledict[wcs.f],
    output:
        read1=config["fq_symlink_dir"] + "/{f}_R1.fastq.gz",
        read2=config["fq_symlink_dir"] + "/{f}_R2.fastq.gz",	
    shell:
        """
        if [ -f {output.read1} ]; then \\rm {output.read1}; fi
        if [ -f {output.read2} ]; then \\rm {output.read2}; fi
        ln -s --relative "{config[raw_fq_dir]}/{params.old_sample_id}_R1.fastq.gz" {output.read1}
        ln -s --relative "{config[raw_fq_dir]}/{params.old_sample_id}_R2.fastq.gz" {output.read2}
        """

rule trimmomatic:
    input:
        read1 = config["fq_symlink_dir"] + "/{fq_id}_R1.fastq.gz",
        read2 = config["fq_symlink_dir"] + "/{fq_id}_R2.fastq.gz",
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
        bam = config["bam_dir"] + "/{fq_id}_raw.bam",
        dedup = temp(config["bam_dir"] + "/{fq_id}_dedup_unsort.bam"),
        sort = config["bam_dir"] + "/{fq_id}_dedup.bam",
        index = config["bam_dir"] + "/{fq_id}_dedup.bam.bai",
    log:
        config["log_dir"] + "/alignment_processing_{fq_id}.log"
    shell:
        """
        sambamba view -t {config[threads]} -S -f bam {input} > {output.bam}
        sambamba markdup -r -t {config[threads]} {output.bam} {output.dedup}
        sambamba sort -t {config[threads]} {output.dedup} -o {output.sort}
        sambamba index -t {config[threads]} {output.sort}
        """

rule fastqc:
    input: 
        raw=config["fq_symlink_dir"] + "/{read_id}_{read}.fastq.gz",
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

rule alignment_qc:
    input:
        config["bam_dir"] + "/{fq_id}_{bam_step}.bam",
    output:
        samstat = config["qc_dir"] + "/{fq_id}_{bam_step}_samstats.txt",
        flagstat = config["qc_dir"] + "/{fq_id}_{bam_step}_flagstat.txt",        
    shell:
        """
        samtools stats {input} > {output.samstat}
        samtools flagstat {input} > {output.flagstat}
        """

rule multiqc:
    input:
        expand(config["qc_dir"] + "/{read_id}_{read}_fastqc.html", read_id = sampledict.keys(), read = ["R1","R2"]),
        expand(config["qc_dir"] + "/{read_id}_proc_{read}_fastqc.html", read_id = sampledict.keys(), read = ["R1","R2"]),
        expand(config["qc_dir"] + "/{read_id}_{bam_step}_samstats.txt", read_id = sampledict.keys(), bam_step = ["raw","dedup"]),
        expand(config["qc_dir"] + "/{read_id}_{bam_step}_flagstat.txt", read_id = sampledict.keys(), bam_step = ["raw","dedup"]),
    params:
        out_dir = config["qc_dir"]
    output:
        config["qc_dir"] + "/all_qc.html"
    shell:
        """
        multiqc {params.out_dir} \
        --force \
        --outdir {params.out_dir} \
        --filename all_qc 
        """

rule downsample_bams:
    input:
        bam = config["bam_dir"] + "/{fq_id}_dedup.bam",
    output:
        config["bam_dir"] + "/{fq_id}_ds{milreads}.bam",
    shell:
        """
        reads=$(echo {wildcards.milreads}000000)
        workflow/scripts/downsample_bam.sh {input} $reads {output}
        """
