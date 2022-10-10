##############################
###   cfDNA WGS Pipeline   ###
##############################

# Read trimming per NCI
rule trimmomatic:
    input:
        read1 = cfdna_wgs_fastq_dir + "/raw/{library_id}_R1.fastq.gz",
        read2 =  cfdna_wgs_fastq_dir + "/raw/{library_id}_R2.fastq.gz",
    params:
        adapter_fasta = config["adapter_fastq"],
	    script = config["cfdna_wgs_script_dir"] + "/trimmomatic_wrapper.sh",
    output:
        read1 = cfdna_wgs_fastq_dir + "/processed/{library_id}_proc_R1.fastq.gz",
        read1_unpr = cfdna_wgs_fastq_dir + "/unpaired/{library_id}_unpr_R1.fastq.gz",
        read2 = cfdna_wgs_fastq_dir + "/processed/{library_id}_proc_R2.fastq.gz",
        read2_unpr = cfdna_wgs_fastq_dir + "/unpaired/{library_id}_unpr_R2.fastq.gz",
    log:
        int = cfdna_wgs_log_dir + "/trimmomatic_trimlog_cfdna_wgs_{library_id}.log",
        main = cfdna_wgs_log_dir + "/trimmomatic_cfdna_wgs_{library_id}.log",
    container:
        config["cfdna_wgs_container"]
    shell:
        """
        {params.script} \
        {input.read1} \
        {input.read2} \
        {params.adapter_fasta} \
        {config[threads]} \
        {output.read1} \
        {output.read1_unpr} \
        {output.read2} \
        {output.read2_unpr} \
        {log.int} \
        &> {log.main}
        """

# FastQC
rule fastqc:
    input:
        raw =  cfdna_wgs_fastq_dir + "/raw/{library_id}_{read}.fastq.gz",
        proc = cfdna_wgs_fastq_dir + "/processed/{library_id}_proc_{read}.fastq.gz",
    params:
        out_dir = cfdna_wgs_qc_dir
    output:
        raw_html = cfdna_wgs_qc_dir + "/{library_id}_{read}_fastqc.html",
        proc_html = cfdna_wgs_qc_dir + "/{library_id}_proc_{read}_fastqc.html",
    log:
        raw = cfdna_wgs_log_dir + "/fastqc_raw_{library_id}_{read}.log",
        proc = cfdna_wgs_log_dir + "/fastqc_proc_{library_id}_{read}.log",
    container:
        config["cfdna_wgs_container"]
    shell:
        """
        fastqc --outdir {params.out_dir} \
        --quiet \
        --threads {config[threads]} {input.raw} &> {log.raw}
        fastqc --outdir {params.out_dir} \
        --quiet \
        --threads {config[threads]} {input.proc} &> {log.proc}
        """

rule index:
    input:
        config["genome_fasta"],
    params:
        out_prefix = genome_ref
    output:
        done = touch(genome_ref)
    container:
        config["cfdna_wgs_container"]
    shell:
        """
        bwa index -p {params.out_prefix} {input}
        """

# BWA alignment
# Post-processing with samblaster and samtools
# Final bam is duplicate marked (NOT removed), location sorted
rule align:
    input:
        ref = genome_ref,
        r1 = cfdna_wgs_fastq_dir + "/processed/{library_id}_proc_R1.fastq.gz",
        r2 = cfdna_wgs_fastq_dir + "/processed/{library_id}_proc_R2.fastq.gz",
    params:
        script = config["cfdna_wgs_script_dir"] + "/align.sh"
    output:
        sort = cfdna_wgs_bam_dir + "/raw/{library_id}.bam",
        index = cfdna_wgs_bam_dir + "/raw/{library_id}.bam.bai",
    log:
        cfdna_wgs_log_dir + "/align_{library_id}.log"
    container:
        config["cfdna_wgs_container"]
    shell:
        """
        {params.script} \
        {input.ref} \
        {input.r1} \
        {input.r2} \
        {config[threads]} \
        {output.sort} &> {log}
	"""

# Alignment samtools QC
rule alignment_qc:
    input:
        cfdna_wgs_bam_dir + "/raw/{library_id}.bam",
    params:
        threads = config["threads"],
    output:
        samstat = cfdna_wgs_qc_dir + "/{library_id}_samstats.txt",
        flagstat = cfdna_wgs_qc_dir + "/{library_id}_flagstat.txt",
    log:
        cfdna_wgs_qc_dir + "/alignment_qc_{library_id}.log",
    container:
        config["cfdna_wgs_container"]
    shell:
        """
        samtools stats -@ {params.threads} {input} > {output.samstat} 2>{log}
        samtools flagstat -@ {params.threads} {input} > {output.flagstat} 2>{log}
        """

# Removes unmapped, not primary, and duplicate reads. Additionally, quality filters by config variable.
rule alignment_filtering:
    input:
        cfdna_wgs_bam_dir + "/raw/{library_id}.bam",
    params:
        script = config["cfdna_wgs_script_dir"] + "/alignment_filtering.sh",
        quality = config["qscore"],
        threads = config["threads"],
    output:
        bam = cfdna_wgs_bam_dir + "/filt/{library_id}_filt.bam",
        bai = cfdna_wgs_bam_dir + "/filt/{library_id}_filt.bam.bai",
    log:
        cfdna_wgs_log_dir + "/{library_id}_alignment_filtering.log",
    container:
        config["cfdna_wgs_container"]
    shell:
        """
        {params.script} \
        {input} \
        {params.quality} \
        {params.threads} \
        {output.bam} &> {log}
        """

# Sequencing depth via Picard
rule picard_collect_wgs_metrics:
    input:
        cfdna_wgs_bam_dir + "/filt/{library_id}_filt.bam",
    params:
        script = config["cfdna_wgs_script_dir"] + "/CollectWgsMetrics_wrapper.sh",
    output:
        cfdna_wgs_qc_dir + "/{library_id}_collect_wgs_metrics.txt",
    log:
        cfdna_wgs_log_dir + "/{library_id}_picard_wgs.log",
    container:
        config["cfdna_wgs_container"]
    shell:
        """
        {config[cfdna_wgs_script_dir]}/CollectWgsMetrics_wrapper.sh \
        {input} \
        {config[picard_jar]} \
        {config[genome_fasta]} \
        {output}
        """

# Fragment sizes by deepTools
rule deeptools_bamprfragmentsize:
    input:
        cfdna_wgs_bam_dir + "/filt/{library_id}_filt.bam",
    params:
        blacklist = config["blacklist"],
        script = config["cfdna_wgs_script_dir"] + "/bamPEFragmentSize_wrapper.sh",
    output:
        cfdna_wgs_qc_dir + "/{library_id}_deeptools_frag_lengths.txt",
    container:
        config["cfdna_wgs_container"]
    shell:
        """
        {params.script} \
        {input} \
        {config[threads]} \
        {params[blacklist]} \
        {output}
        """

rule cfdna_wgs_multiqc:
    input:
        expand(cfdna_wgs_qc_dir + "/{library_id}_{read}_fastqc.html", library_id = LIBRARIES, read = ["R1","R2"]),
        expand(cfdna_wgs_qc_dir + "/{library_id}_proc_{read}_fastqc.html", library_id = LIBRARIES, read = ["R1","R2"]),
        expand(cfdna_wgs_qc_dir + "/{library_id}_samstats.txt", library_id = LIBRARIES),
        expand(cfdna_wgs_qc_dir + "/{library_id}_flagstat.txt", library_id = LIBRARIES),
        expand(cfdna_wgs_qc_dir + "/{library_id}_deeptools_frag_lengths.txt", library_id = LIBRARIES),
        expand(cfdna_wgs_qc_dir + "/{library_id}_deeptools_frag_lengths.txt", library_id = LIBRARIES),
        expand(cfdna_wgs_qc_dir + "/{library_id}_collect_wgs_metrics.txt", library_id = LIBRARIES),
    params:
        out_dir = cfdna_wgs_qc_dir
    output:
        cfdna_wgs_qc_dir + "/all_qc_data/multiqc_fastqc.txt",
        cfdna_wgs_qc_dir + "/all_qc_data/multiqc_samtools_stats.txt",
        cfdna_wgs_qc_dir + "/all_qc_data/multiqc_samtools_flagstat.txt",
	    cfdna_wgs_qc_dir + "/all_qc_data/multiqc_picard_wgsmetrics.txt",
    container:
        config["cfdna_wgs_container"]
    shell:
        """
        multiqc {params.out_dir} \
        --force \
        --outdir {params.out_dir} \
        --filename all_qc
        """

rule aggregate_frag:
    input:
        expand(cfdna_wgs_qc_dir + "/{library_id}_deeptools_frag_lengths.txt", library_id = LIBRARIES),
    params:
        script = config["cfdna_wgs_script_dir"] + "/aggregate_frag.sh",
    output:
        cfdna_wgs_qc_dir + "/all_frag.tsv",
    log:
        cfdna_wgs_log_dir + "/aggregate_frag.err",
    container:
        config["cfdna_wgs_container"]
    shell:
        """
        awk 'FNR>2' {input} > {output} 2> {log}
        """

#  Notes:
#  This makes an aggregate table of QC values. The subsequent downsampling
#  step only runs if read numbers are above a certain threshold. See also
#  the int_test.smk for function using this output table.
#

checkpoint make_qc_tbl:
    input:
        fq = cfdna_wgs_qc_dir + "/all_qc_data/multiqc_fastqc.txt",
        sam = cfdna_wgs_qc_dir + "/all_qc_data/multiqc_samtools_stats.txt",
        flag = cfdna_wgs_qc_dir + "/all_qc_data/multiqc_samtools_flagstat.txt",
	    picard = cfdna_wgs_qc_dir + "/all_qc_data/multiqc_picard_wgsmetrics.txt",
        deeptools = cfdna_wgs_qc_dir + "/all_frag.tsv",
    params:
        script = config["cfdna_wgs_script_dir"] + "/make_qc_tbl.R"
    output:
        cfdna_wgs_qc_dir + "/read_qc.tsv",
    log:
        cfdna_wgs_log_dir + "/read_qc.log"
    container:
        config["cfdna_wgs_container"]
    shell:
        """
        Rscript {params.script} \
        {input.fq} \
        {input.sam} \
        {input.flag} \
        {input.picard} \
        {input.deeptools} \
        {output} \
        >& {log}
        """

# Alignment downsampling
#  Note: Used for all rule input "get_ds_candidates". See that function in
#  workflow/int_test.smk

rule downsample_bams:
    input:
        cfdna_wgs_bam_dir + "/filt/{library_id}_filt.bam",
    output:
        cfdna_wgs_bam_dir + "/ds/{library_id}_ds{milreads}.bam",
    log:
        cfdna_wgs_log_dir + "/downsample_bam_{library_id}_{milreads}.err"
    container:
        config["cfdna_wgs_container"]
    shell:
        """
        {config[cfdna_wgs_script_dir]}/downsample_bam.sh {input} {wildcards.milreads} {output} {config[threads]} 2>{log}
        """
