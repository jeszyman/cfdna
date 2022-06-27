# Read trimming per NCI
rule trimmomatic:
    input:
        read1 = config["data_dir"] + "/fastq/raw/{library_id}_R1.fastq.gz",
        read2 = config["data_dir"] + "/fastq/raw/{library_id}_R2.fastq.gz",
    params:
        adapter_fasta = config["adapter_fastq"],
	script = config["cfdna_wgs_script_dir"] + "/trimmomatic_wrapper.sh",
    output:
        read1 = config["data_dir"] + "/fastq/processed/{library_id}_proc_R1.fastq.gz",
        read1_unpr = config["data_dir"] + "/fastq/unpaired/{library_id}_unpr_R1.fastq.gz",
        read2 = config["data_dir"] + "/fastq/processed/{library_id}_proc_R2.fastq.gz",
        read2_unpr = config["data_dir"] + "/fastq/unpaired/{library_id}_unpr_R2.fastq.gz",
    log:
        int = config["data_dir"] + "/logs/trimmomatic_trimlog_{library_id}.log",
        main = config["data_dir"] + "/logs/trimmomatic_{library_id}.log",
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
        raw =  config["data_dir"] + "/fastq/raw/{library_id}_{read}.fastq.gz",
        proc = config["data_dir"] + "/fastq/processed/{library_id}_proc_{read}.fastq.gz",
    params:
        out_dir = config["data_dir"] + "/qc",
    output:
        raw_html = config["data_dir"] + "/qc/{library_id}_{read}_fastqc.html",
        proc_html = config["data_dir"] + "/qc/{library_id}_proc_{read}_fastqc.html",
    log:
        raw = config["data_dir"] + "/logs/fastqc_raw_{library_id}_{read}.log",
        proc = config["data_dir"] + "/logs/fastqc_proc_{library_id}_{read}.log",
    shell:
        """
        fastqc --outdir {params.out_dir} \
        --quiet \
        --threads {config[threads]} {input.raw} &> {log}
        fastqc --outdir {params.out_dir} \
        --quiet \
        --threads {config[threads]} {input.proc} &> {log}
        """

rule index:
    input:
        config["genome_fasta"],
    params:
        out_prefix = genome_ref
    output:
        done = touch(genome_ref)
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
        r1 = config["data_dir"] + "/fastq/processed/{library_id}_proc_R1.fastq.gz",
        r2 = config["data_dir"] + "/fastq/processed/{library_id}_proc_R2.fastq.gz",
    params:
        script = config["cfdna_wgs_script_dir"] + "/align.sh"
    output:
        sort = config["data_dir"] + "/bam/raw/{library_id}.bam",
        index = config["data_dir"] + "/bam/raw/{library_id}.bam.bai",
    log:
        config["data_dir"] + "/logs/align_{library_id}.log"
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
        config["data_dir"] + "/bam/raw/{library_id}.bam",
    params:
        threads = config["threads"],
    output:
        samstat = config["data_dir"] + "/qc/{library_id}_samstats.txt",
        flagstat = config["data_dir"] + "/qc/{library_id}_flagstat.txt",
    log:
        config["data_dir"] + "/logs/alignment_qc_{library_id}.log",
    shell:
        """
        samtools stats -@ {params.threads} {input} > {output.samstat} 2>{log}
        samtools flagstat -@ {params.threads} {input} > {output.flagstat} 2>{log}
        """

# Removes unmapped, not primary, and duplicate reads. Additionally, quality filters by config variable.
rule alignment_filtering:
    input:
        config["data_dir"] + "/bam/raw/{library_id}.bam",
    params:
        script = config["cfdna_wgs_script_dir"] + "/alignment_filtering.sh",
        quality = config["qscore"],
        threads = config["threads"],
    output:
        bam = config["data_dir"] + "/bam/filt/{library_id}_filt.bam",
        bai = config["data_dir"] + "/bam/filt/{library_id}_filt.bam.bai",
    log:
        config["data_dir"] + "/logs/{library_id}_alignment_filtering.log",
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
        config["data_dir"] + "/bam/filt/{library_id}_filt.bam",
    params:
        script = config["cfdna_wgs_script_dir"] + "/CollectWgsMetrics_wrapper.sh",
    output:
        config["data_dir"] + "/qc/{library_id}_collect_wgs_metrics.txt",
    log:
        config["data_dir"] + "/logs/{library_id}_picard_wgs.log",
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
        config["data_dir"] + "/bam/filt/{library_id}_filt.bam",
    params:
        blacklist = config["blacklist"],
        script = config["cfdna_wgs_script_dir"] + "/bamPEFragmentSize_wrapper.sh",
    output:
        config["data_dir"] + "/qc/{library_id}_deeptools_frag_lengths.txt",
    shell:
        """
        {params.script} \
        {input} \
        {config[threads]} \
        {params[blacklist]} \
        {output}
        """

rule multiqc:
    input:
        expand(config["data_dir"] + "/qc/{library_id}_{read}_fastqc.html", library_id = LIBRARY_IDS, read = ["R1","R2"]),
        expand(config["data_dir"] + "/qc/{library_id}_proc_{read}_fastqc.html", library_id = LIBRARY_IDS, read = ["R1","R2"]),
        expand(config["data_dir"] + "/qc/{library_id}_samstats.txt", library_id = LIBRARY_IDS),
        expand(config["data_dir"] + "/qc/{library_id}_flagstat.txt", library_id = LIBRARY_IDS),
        expand(config["data_dir"] + "/qc/{library_id}_deeptools_frag_lengths.txt", library_id = LIBRARY_IDS),
        expand(config["data_dir"] + "/qc/{library_id}_deeptools_frag_lengths.txt", library_id = LIBRARY_IDS),
        expand(config["data_dir"] + "/qc/{library_id}_collect_wgs_metrics.txt", library_id = LIBRARY_IDS),
    params:
        out_dir = config["data_dir"] + "/qc"
    output:
        config["data_dir"] + "/qc/all_qc_data/multiqc_fastqc.txt",
        config["data_dir"] + "/qc/all_qc_data/multiqc_samtools_stats.txt",
        config["data_dir"] + "/qc/all_qc_data/multiqc_samtools_flagstat.txt",
	config["data_dir"] + "/qc/all_qc_data/multiqc_picard_wgsmetrics.txt",
    shell:
        """
        multiqc {params.out_dir} \
        --force \
        --outdir {params.out_dir} \
        --filename all_qc
        """

rule aggregate_frag:
    input:
        expand(config["data_dir"] + "/qc/{library_id}_deeptools_frag_lengths.txt", library_id = LIBRARY_IDS),
    params:
        script = config["cfdna_wgs_script_dir"] + "/aggregate_frag.sh",
    output:
        config["data_dir"] + "/qc/all_frag.tsv",
    log:
        config["data_dir"] + "/logs/aggregate_frag.err",
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
        fq = config["data_dir"] + "/qc/all_qc_data/multiqc_fastqc.txt",
        sam = config["data_dir"] + "/qc/all_qc_data/multiqc_samtools_stats.txt",
        flag = config["data_dir"] + "/qc/all_qc_data/multiqc_samtools_flagstat.txt",
	picard = config["data_dir"] + "/qc/all_qc_data/multiqc_picard_wgsmetrics.txt",
        deeptools = config["data_dir"] + "/qc/all_frag.tsv",
    params:
        script = config["cfdna_wgs_script_dir"] + "/make_qc_tbl.R"
    output:
        config["data_dir"] + "/qc/read_qc.tsv",
    log:
        config["data_dir"] + "/logs/read_qc.log"
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
        config["data_dir"] + "/bam/filt/{library_id}_filt.bam",
    output:
        config["data_dir"] + "/bam/ds/{library_id}_ds{milreads}.bam",
    log:
        config["data_dir"] + "/logs/downsample_bam_{library_id}_{milreads}.err"
    shell:
        """
        {config[cfdna_wgs_script_dir]}/downsample_bam.sh {input} {wildcards.milreads} {output} 2>{log}
        """
