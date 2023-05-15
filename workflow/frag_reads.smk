##################################################################
###   Integration testing snakefile for WGS cfDNA Processing   ###
##################################################################

rule frag_index:
    benchmark: "{bench_dir}/frag_index.benchmark.txt",
    input: genome_fasta,
    log: "{log_dir}/frag_index.log",
    output: "{data_dir}/ref/{fasta_base}.sa",
    params:
        out_prefix = "{bwa_dir}/{fasta_base}",
        script = "{frag_script_dir}/frag_index.sh",
    shell:
        """
        {params.script} {input} {params.out_prefix} &> {log}
        """

# Adapter-trim and QC reads with fastp
rule frag_fastp:
    benchmark: "{bench_dir}/{{library}}_frag_fastp.benchmark.txt",
    input:
        read1 = "{frag_fastq_dir/{{library}}_raw_R1.fastq.gz",
        read2 = "{frag_fastq_dir/{{library}}_raw_R2.fastq.gz",
    log:
        cmd = "{log_dir/{{library}}_frag_fastp.log",
        html = "{log_dir/{{library}}_frag_fastp.html",
        json = "{log_dir/{{library}}_frag_fastp.json",
    output:
        read1 = "{frag_fastq_dir}/{{library}}_processed_R1.fastq.gz",
        read2 = "{frag_fastq_dir}/{{library}}_processed_R2.fastq.gz",
        failed = "{frag_fastq_dir}/{{library}}_failed_fastp.fastq.gz",
        unpaired1 = "{frag_fastq_dir}/{{library}}_unpaired_R1.fastq.gz",
        unpaired2 = "{frag_fastq_dir}/{{library}}_unpaired_R2.fastq.gz",
    params:
        script = "{frag_script_dir}/fastp.sh",
        threads = frag_threads,
    resources:
        mem_mb = 500,
    shell:
        """
        {params.script} \
        {input.read1} \
        {input.read2} \
        {log.html} \
        {log.json} \
        {output.read1} \
        {output.read2} \
        {output.failed} \
        {output.unpaired1} \
        {output.unpaired2} \
        {params.threads} &> {log.cmd}
        """

# Align reads with BWA
rule frag_align:
    benchmark: benchdir + "/{library}_frag_align.benchmark.txt",
    input:
        ref = "{data_dir}/ref/{fasta_base}",
        read1 = frag_fastqs + "/{library}_processed_R1.fastq.gz",
        read2 = frag_fastqs + "/{library}_processed_R2.fastq.gz",
    log: logdir + "/{library}_frag_align.log",
    output:
        sort = frag_bams + "/{library}_raw.bam",
        index = frag_bams + "/{library}_raw.bam.bai",
    params:
        script = "{frag_script_dir}/align.sh",
        threads = 4,
    resources:
        mem_mb = 500,
    shell:
        """
        {params.script} \
        {input.ref} \
        {input.read1} \
        {input.read2} \
        {params.threads} \
        {output.sort} &> {log}
        """

# Remove PCR duplicates from aligned reads
rule frag_dedup:
    benchmark: benchdir + "/{library}_frag_dedup.benchmark.txt",
    input: frag_bams + "/{library}_raw.bam",
    log: logdir + "/{library}_frag_dedup.log",
    output: frag_bams + "/{library}_dedup.bam",
    params:
        script = "{frag_script_dir}/dedup.sh",
        threads = frag_threads,
    shell:
        """
        {params.script} \
        {input} \
        {output} \
        {params.threads} &> {log}
        """

# Filter de-duplicated alignments.
# Remove unmapped, not primary, and duplicate reads. Additional location filter by config bedfile variable.

checkpoint frag_filter_alignment:
    benchmark: benchdir + "/{library}_frag_filter_alignment.benchmark.txt",
    input: frag_bams + "/{library}_dedup.bam",
    log: logdir + "/{library}_frag_filter_alignment.log",
    output: frag_bams + "/{library}_filt.bam",
    params:
        script = "{frag_script_dir}/filter_alignment.sh",
        threads = frag_threads,
    shell:
        """
        {params.script} \
        {input} \
        {params.threads} \
        {output} &> {log}
        """

# Get read quality by FASTQC
rule frag_fastqc:
    benchmark: benchdir+ "/{library}_{processing}_{read}_frag_fastqc.benchmark.txt",
    input: frag_fastqs + "/{library}_{processing}_{read}.fastq.gz",
    log: logdir + "/{library}_{processing}_{read}_frag_fastqc.log",
    output:
        qcdir + "/{library}_{processing}_{read}_fastqc.html",
        qcdir + "/{library}_{processing}_{read}_fastqc.zip",
    params:
        outdir = qcdir,
        script = "{frag_script_dir}/fastqc.sh",
        threads = frag_threads,
    shell:
        """
        {params.script} \
        {input} \
        {params.outdir} \
        {params.threads} &> {log}
        """

# Get alignment QC using samtools
rule frag_alignment_qc:
    input: frag_bams + "/{library}_{processing}.bam",
    log:
        flagstat = logdir + "/{library}_{processing}_flagstat_frag_alignment_qc.log",
        samstat = logdir + "/{library}_{processing}_samstats_frag_alignment_qc.log",
    output:
        flagstat = qcdir + "/{library}_{processing}_flagstat.txt",
        samstat = qcdir + "/{library}_{processing}_samstats.txt",
    params:
        script = "{frag_script_dir}/alignment_qc.sh",
        threads = frag_threads,
    shell:
        """
        {params.script} \
        {input} \
        {log.flagstat} \
        {log.samstat} \
        {output.flagstat} \
        {output.samstat} \
        {params.threads}
        """

# Sequencing depth metrics via Picard
rule frag_picard_depth:
    benchmark: benchdir + "/{library}_frag_picard_depth.benchmark.txt",
    input: frag_bams + "/{library}_filt.bam",
    log: logdir + "/{library}_frag_picard_depth.log",
    output: qcdir + "/{library}_picard_depth.txt",
    params:
        script = "{frag_script_dir}/picard_depth.sh",
        threads = frag_threads,
    shell:
        """
        {params.script} \
        {input} \
        {config[picard_jar]} \
        {config[genome_fasta]} \
        {output}
        """

# Get fragment sizes using deepTools
rule frag_bampefragsize:
    benchmark: benchdir + "/frag_bampefragsize.benchmark.txt",
    input: expand(frag_bams + "/{library}_filt.bam", library = FRAG_LIBS),
    log: logdir + "/frag_bampefragsize.log",
    output:
        raw = qcdir + "/deeptools_frag_lengths.txt",
        hist = qcdir + "/deeptools_frag_lengths.png",
    params:
        blacklist = config["blklist"],
        script = "{frag_script_dir}/bampefragsize.sh",
        threads = frag_threads,
    shell:
        """
        {params.script} \
        "{input}" \
        {log} \
        {output.hist} \
        {output.raw} \
        {params.blacklist} \
        {params.threads}
        """

# Make deeptools bamCoverage bedfile
rule frag_bamcoverage:
    benchmark: benchdir + "/{library}_frag_bamcoverage.benchmark.txt",
    input: frag_bams + "/{library}_filt.bam",
    log: logdir + "/{library}_frag_bamcoverage.log",
    output: qcdir + "/{library}_bamcoverage.bg",
    params:
        bin = "10000",
        blacklist = config["blklist"],
        script = "{frag_script_dir}/bamcoverage.sh",
        threads = frag_threads,
    shell:
        """
        {params.script} \
        {input} \
        {output} \
        {params.bin} \
        {params.blacklist} \
        {params.threads} &> {log}
        """

# Make deepTools plotCoverage coverage maps for all filtered bams
rule frag_plotcoverage:
    benchmark: benchdir + "/frag_plotcoverage.benchmark.txt",
    input: expand(frag_bams + "/{library}_filt.bam", library = FRAG_LIBS),
    log: logdir + "/frag_plotcoverage.log",
    output:
        raw = qcdir + "/frag_coverage.tsv",
        plot = qcdir + "/frag_coverage.pd",
    params:
        blacklist = config["blklist"],
        script = "{frag_script_dir}/plotcoverage.sh",
        threads = frag_threads,
    shell:
        """
        {params.script} \
        "{input}" \
        {params.blacklist} \
        {params.threads} \
        {output.raw} \
        {output.plot} &> {log}
        """

# Aggregate QC files using MultiQC
rule frag_multiqc:
    benchmark: benchdir + "/frag_multiqc.benchmark.txt",
    input:
        expand(logdir + "/{library}_frag_fastp.json", library = FRAG_LIBS),
        expand(qcdir + "/{library}_{processing}_{read}_fastqc.zip", library = FRAG_LIBS, processing = ["raw", "processed", "unpaired"], read = ["R1","R2"]),
        expand(qcdir + "/{library}_{processing}_samstats.txt", library = FRAG_LIBS, processing = ["raw","filt"]),
        expand(qcdir + "/{library}_{processing}_flagstat.txt", library = FRAG_LIBS, processing = ["raw","filt"]),
        expand(qcdir + "/{library}_picard_depth.txt", library = FRAG_LIBS),
        qcdir + "/deeptools_frag_lengths.txt",
        qcdir + "/frag_coverage.tsv",
    log: logdir + "/frag_multiqc.log",
    output:
        qcdir + "/frag_multiqc.html",
        qcdir + "/frag_multiqc_data/multiqc_fastqc.txt",
        qcdir + "/frag_multiqc_data/multiqc_samtools_stats.txt",
        qcdir + "/frag_multiqc_data/multiqc_picard_wgsmetrics.txt",
        qcdir + "/frag_multiqc_data/multiqc_samtools_flagstat.txt",
    params:
        out_dir = qcdir,
        out_name = "frag_multiqc",
        script = "{frag_script_dir}/multiqc.sh",
        threads = frag_threads,
    shell:
        """
        {params.script} \
        "{input}" \
        {params.out_name} \
        {params.out_dir} &> {log}
        """

# Make a tab-separated aggregate QC table
checkpoint frag_make_qc_tsv:
    benchmark: benchdir + "/frag_make_qc_tsv.benchmark.txt",
    input:
        fq = qcdir + "/frag_multiqc_data/multiqc_fastqc.txt",
        mqsam = qcdir + "/frag_multiqc_data/multiqc_samtools_stats.txt",
        mqflag = qcdir + "/frag_multiqc_data/multiqc_samtools_flagstat.txt",
        picard = qcdir + "/frag_multiqc_data/multiqc_picard_wgsmetrics.txt",
        deeptools_frag = qcdir + "/deeptools_frag_lengths.txt",
        deeptools_cov = qcdir + "/frag_coverage.tsv",
    log: logdir + "/frag_make_qc_tsv.log",
    output:
        readqc = qcdir + "/frag_read_qc.tsv",
        fraglen = qcdir + "/frag_frag_len.tsv",
    params:
        script = "{frag_script_dir}/make_qc_tsv.R",
    shell:
        """
        Rscript {params.script} \
        {input.fq} \
        {input.mqsam} \
        {input.mqflag} \
        {input.picard} \
        {input.deeptools_frag} \
        {input.deeptools_cov} \
        {output.readqc} \
        {output.fraglen} >& {log}
        """
