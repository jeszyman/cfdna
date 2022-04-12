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
        expand(config["qc_dir"] + "/{read_id}_{read}_fastqc.html", read_id = sampledict.keys(), read = ["R1","R2"]),
        expand(config["qc_dir"] + "/{read_id}_proc_{read}_fastqc.html", read_id = sampledict.keys(), read = ["R1","R2"]),
	expand(config["bam_dir"] + "/{read_id}_dedup_sort.bam", read_id = sampledict.keys())
	
checkpoint rename:
    params:
        old_sample_id=lambda wcs: sampledict[wcs.f],
        read1=config["fq_symlink_dir"] + "/{f}_R1.fastq.gz",
        read2=config["fq_symlink_dir"] + "/{f}_R2.fastq.gz",
    output:
        config["fq_symlink_dir"]
    shell:
        """
        if [ -f {params.read1} ]; then \\rm {params.read1}; fi
        if [ -f {params.read2} ]; then \\rm {params.read2}; fi
        ln -s --relative "{config[raw_fq_dir]}/{params.old_sample_id}_R1.fastq.gz" {params.read1}
        ln -s --relative "{config[raw_fq_dir]}/{params.old_sample_id}_R2.fastq.gz" {params.read2}
        """

def symlink_out(wildcards):
    checkpoint_output = checkpoints.get(**wildcards).output[0]    
    file_names = expand("outputs/cd-hit95/{mag}.cdhit95.faa", 
                        mag = glob_wildcards(os.path.join(checkpoint_output, "{mag}.fa.cdbg_ids.reads.hardtrim.fa.gz.plass.cdhit.fa.clean.cut.dup")).mag)
    return file_names

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
