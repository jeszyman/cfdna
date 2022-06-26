# Read trimming per NCI
rule trimmomatic:
    input:
        read1 = config["data_dir"] + "/fastq/raw/{library_id}_R1.fastq.gz",
        read2 = config["data_dir"] + "/fastq/raw/{library_id}_R2.fastq.gz",
    params:
        adapter_fasta = config["adapter_fastq"],
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
        {config[cfdna_wgs_script_dir]}/trimmomatic_wrapper.sh \
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
rule align:
    input:
        bwa_index_done = genome_ref,
        read1 = config["data_dir"] + "/fastq/processed/{library_id}_proc_R1.fastq.gz",
        read2 = config["data_dir"] + "/fastq/processed/{library_id}_proc_R2.fastq.gz",
    output:
        config["data_dir"] + "/bam/{library_id}.sam",
    log:
        config["data_dir"] + "/logs/align_{library_id}.log"
    shell:
        """
        bwa mem -M -t 4 {input.bwa_index_done} {input.read1} {input.read2} > {output}
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

# Alignment deduplication and sorting
rule alignment_processing:
    input:
        config["data_dir"] + "/bam/{library_id}.sam",
    output:
        bam = config["data_dir"] + "/bam/{library_id}_raw.bam",
        dedup = temp(config["data_dir"] + "/bam/{library_id}_dedup_unsort.bam"),
        sort = config["data_dir"] + "/bam/{library_id}_dedup.bam",
        index = config["data_dir"] + "/bam/{library_id}_dedup.bam.bai",
    log:
        config["data_dir"] + "/logs/alignment_processing_{library_id}.log"
    shell:
        """
        {config[cfdna_wgs_script_dir]}/alignment_processing.sh \
        {input} \
        {config[threads]} \
        {output.bam} \
        {output.dedup} \
        {output.sort} \
        {output.index} \
        &> {log}
        """

# Alignment samtools QC
rule alignment_qc:
    input:
        config["data_dir"] + "/bam/{library_id}_{bam_step}.bam",
    output:
        samstat = config["data_dir"] + "/qc/{library_id}_{bam_step}_samstats.txt",
        flagstat = config["data_dir"] + "/qc/{library_id}_{bam_step}_flagstat.txt",
    log:
        config["data_dir"] + "/logs/alignment_qc_{library_id}_{bam_step}.err",
    shell:
        """
        samtools stats {input} > {output.samstat} 2>{log}
        samtools flagstat {input} > {output.flagstat} 2>{log}
        """

# Sequencing depth via Picard
rule picard_collect_wgs_metrics:
    input:
        config["data_dir"] + "/bam/{library_id}_dedup.bam",
    output:
        config["data_dir"] + "/qc/{library_id}_collect_wgs_metrics.txt",
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
        config["data_dir"] + "/bam/{library_id}_dedup.bam",
    params:
        blacklist = config["blacklist"],
    output:
        config["data_dir"] + "/qc/{library_id}_deeptools_frag_lengths.txt",
    shell:
        """
        {config[cfdna_wgs_script_dir]}/bamPEFragmentSize_wrapper.sh \
        {input} \
        {config[threads]} \
        {params[blacklist]} \
        {output}
        """

rule multiqc:
    input:
        expand(config["data_dir"] + "/qc/{library_id}_{read}_fastqc.html", library_id = LIBRARY_IDS, read = ["R1","R2"]),
        expand(config["data_dir"] + "/qc/{library_id}_proc_{read}_fastqc.html", library_id = LIBRARY_IDS, read = ["R1","R2"]),
        expand(config["data_dir"] + "/qc/{library_id}_{bam_step}_samstats.txt", library_id = LIBRARY_IDS, bam_step= ["dedup","raw"]),
        expand(config["data_dir"] + "/qc/{library_id}_{bam_step}_flagstat.txt", library_id = LIBRARY_IDS, bam_step =["dedup","raw"]),
    params:
        out_dir = config["data_dir"] + "/qc"
    output:
        config["data_dir"] + "/qc/all_qc.html",
        config["data_dir"] + "/qc/all_qc_data/multiqc_samtools_stats.txt",
    shell:
        """
        multiqc {params.out_dir} \
        --force \
        --outdir {params.out_dir} \
        --filename all_qc
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
        {output} \
        >& {log}
        """

# Alignment downsampling
rule downsample_bams:
    input:
        config["data_dir"] + "/bam/{library_id}_dedup.bam",
    output:
        config["data_dir"] + "/bam/{library_id}_ds{milreads}.bam",
    log:
        config["data_dir"] + "/logs/downsample_bam_{library_id}_{milreads}.err"
    shell:
        """
        {config[cfdna_wgs_script_dir]}/downsample_bam.sh {input} {wildcards.milreads}000000 {output} 2>{log}
        """
