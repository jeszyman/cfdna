#########1#########2#########3#########4#########5#########6#########7#########8
#                                                                              #
#                    Basic Read and Alignment Processing of                    #
#                    Cell-free DNA Whole Genome Sequencing                     #
#                                                                              #
#########1#########2#########3#########4#########5#########6#########7#########8

rule frag_index:
    benchmark: benchdir + "/frag_index.benchmark.txt",
    input: genome_fasta,
    log: logdir + "/frag_index.log",
    output: f"{data_dir}/ref/{fasta_base}.sa",
    params:
        out_prefix = f"{bwa_dir}/{fasta_base}",
        script = f"{frag_script_dir}/frag_index.sh",
    shell:
        """
        {params.script} {input} {params.out_prefix} &> {log}
        """

# Adapter-trim and QC reads with fastp
rule cfdna_wgs_fastp:
    benchmark: benchdir + "/{library}_cfdna_wgs_fastp.benchmark.txt",
    input:
        read1 = cfdna_wgs_fastqs + "/{library}_raw_R1.fastq.gz",
        read2 = cfdna_wgs_fastqs + "/{library}_raw_R2.fastq.gz",
    log:
        cmd = logdir + "/{library}_cfdna_wgs_fastp.log",
        html = logdir + "/{library}_cfdna_wgs_fastp.html",
        json = logdir + "/{library}_cfdna_wgs_fastp.json",
    output:
        read1 = cfdna_wgs_fastqs + "/{library}_processed_R1.fastq.gz",
        read2 = cfdna_wgs_fastqs + "/{library}_processed_R2.fastq.gz",
        failed = cfdna_wgs_fastqs + "/{library}_failed_fastp.fastq.gz",
        unpaired1 = cfdna_wgs_fastqs + "/{library}_unpaired_R1.fastq.gz",
        unpaired2 = cfdna_wgs_fastqs + "/{library}_unpaired_R2.fastq.gz",
    params:
        script = cfdna_wgs_scriptdir + "/fastp.sh",
        threads = cfdna_wgs_threads,
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
rule cfdna_wgs_align:
    benchmark: benchdir + "/{library}_cfdna_wgs_align.benchmark.txt",
    input:
        ref = "{data_dir}/ref/{fasta_base}",
        read1 = cfdna_wgs_fastqs + "/{library}_processed_R1.fastq.gz",
        read2 = cfdna_wgs_fastqs + "/{library}_processed_R2.fastq.gz",
    log: logdir + "/{library}_cfdna_wgs_align.log",
    output:
        sort = cfdna_wgs_bams + "/{library}_raw.bam",
        index = cfdna_wgs_bams + "/{library}_raw.bam.bai",
    params:
        script = cfdna_wgs_scriptdir + "/align.sh",
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
rule cfdna_wgs_dedup:
    benchmark: benchdir + "/{library}_cfdna_wgs_dedup.benchmark.txt",
    input: cfdna_wgs_bams + "/{library}_raw.bam",
    log: logdir + "/{library}_cfdna_wgs_dedup.log",
    output: cfdna_wgs_bams + "/{library}_dedup.bam",
    params:
        script = cfdna_wgs_scriptdir + "/dedup.sh",
        threads = cfdna_wgs_threads,
    shell:
        """
        {params.script} \
        {input} \
        {output} \
        {params.threads} &> {log}
        """

# Filter de-duplicated alignments.
# Remove unmapped, not primary, and duplicate reads. Additional location filter by config bedfile variable.

checkpoint cfdna_wgs_filter_alignment:
    benchmark: benchdir + "/{library}_cfdna_wgs_filter_alignment.benchmark.txt",
    input: cfdna_wgs_bams + "/{library}_dedup.bam",
    log: logdir + "/{library}_cfdna_wgs_filter_alignment.log",
    output: cfdna_wgs_bams + "/{library}_filt.bam",
    params:
        script = cfdna_wgs_scriptdir + "/filter_alignment.sh",
        threads = cfdna_wgs_threads,
    shell:
        """
        {params.script} \
        {input} \
        {params.threads} \
        {output} &> {log}
        """

# Get read quality by FASTQC
rule cfdna_wgs_fastqc:
    benchmark: benchdir+ "/{library}_{processing}_{read}_cfdna_wgs_fastqc.benchmark.txt",
    input: cfdna_wgs_fastqs + "/{library}_{processing}_{read}.fastq.gz",
    log: logdir + "/{library}_{processing}_{read}_cfdna_wgs_fastqc.log",
    output:
        qcdir + "/{library}_{processing}_{read}_fastqc.html",
        qcdir + "/{library}_{processing}_{read}_fastqc.zip",
    params:
        outdir = qcdir,
        script = cfdna_wgs_scriptdir + "/fastqc.sh",
        threads = cfdna_wgs_threads,
    shell:
        """
        {params.script} \
        {input} \
        {params.outdir} \
        {params.threads} &> {log}
        """

# Get alignment QC using samtools
rule cfdna_wgs_alignment_qc:
    input: cfdna_wgs_bams + "/{library}_{processing}.bam",
    log:
        flagstat = logdir + "/{library}_{processing}_flagstat_cfdna_wgs_alignment_qc.log",
        samstat = logdir + "/{library}_{processing}_samstats_cfdna_wgs_alignment_qc.log",
    output:
        flagstat = qcdir + "/{library}_{processing}_flagstat.txt",
        samstat = qcdir + "/{library}_{processing}_samstats.txt",
    params:
        script = cfdna_wgs_scriptdir + "/alignment_qc.sh",
        threads = cfdna_wgs_threads,
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
rule cfdna_wgs_picard_depth:
    benchmark: benchdir + "/{library}_cfdna_wgs_picard_depth.benchmark.txt",
    input: cfdna_wgs_bams + "/{library}_filt.bam",
    log: logdir + "/{library}_cfdna_wgs_picard_depth.log",
    output: qcdir + "/{library}_picard_depth.txt",
    params:
        script = cfdna_wgs_scriptdir + "/picard_depth.sh",
        threads = cfdna_wgs_threads,
    shell:
        """
        {params.script} \
        {input} \
        {config[picard_jar]} \
        {config[genome_fasta]} \
        {output}
        """

# Get fragment sizes using deepTools
rule cfdna_wgs_bampefragsize:
    benchmark: benchdir + "/cfdna_wgs_bampefragsize.benchmark.txt",
    input: expand(cfdna_wgs_bams + "/{library}_filt.bam", library = FRAG_LIBS),
    log: logdir + "/cfdna_wgs_bampefragsize.log",
    output:
        raw = qcdir + "/deeptools_frag_lengths.txt",
        hist = qcdir + "/deeptools_frag_lengths.png",
    params:
        blacklist = config["blklist"],
        script = cfdna_wgs_scriptdir + "/bampefragsize.sh",
        threads = cfdna_wgs_threads,
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
rule cfdna_wgs_bamcoverage:
    benchmark: benchdir + "/{library}_cfdna_wgs_bamcoverage.benchmark.txt",
    input: cfdna_wgs_bams + "/{library}_filt.bam",
    log: logdir + "/{library}_cfdna_wgs_bamcoverage.log",
    output: qcdir + "/{library}_bamcoverage.bg",
    params:
        bin = "10000",
        blacklist = config["blklist"],
        script = cfdna_wgs_scriptdir + "/bamcoverage.sh",
        threads = cfdna_wgs_threads,
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
rule cfdna_wgs_plotcoverage:
    benchmark: benchdir + "/cfdna_wgs_plotcoverage.benchmark.txt",
    input: expand(cfdna_wgs_bams + "/{library}_filt.bam", library = FRAG_LIBS),
    log: logdir + "/cfdna_wgs_plotcoverage.log",
    output:
        raw = qcdir + "/cfdna_wgs_coverage.tsv",
        plot = qcdir + "/cfdna_wgs_coverage.pdf",
    params:
        blacklist = config["blklist"],
        script = cfdna_wgs_scriptdir + "/plotcoverage.sh",
        threads = cfdna_wgs_threads,
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
rule cfdna_wgs_multiqc:
    benchmark: benchdir + "/cfdna_wgs_multiqc.benchmark.txt",
    input:
        expand(logdir + "/{library}_cfdna_wgs_fastp.json", library = FRAG_LIBS),
        expand(qcdir + "/{library}_{processing}_{read}_fastqc.zip", library = FRAG_LIBS, processing = ["raw", "processed", "unpaired"], read = ["R1","R2"]),
        expand(qcdir + "/{library}_{processing}_samstats.txt", library = FRAG_LIBS, processing = ["raw","filt"]),
        expand(qcdir + "/{library}_{processing}_flagstat.txt", library = FRAG_LIBS, processing = ["raw","filt"]),
        expand(qcdir + "/{library}_picard_depth.txt", library = FRAG_LIBS),
        qcdir + "/deeptools_frag_lengths.txt",
        qcdir + "/cfdna_wgs_coverage.tsv",
    log: logdir + "/cfdna_wgs_multiqc.log",
    output:
        qcdir + "/cfdna_wgs_multiqc.html",
        qcdir + "/cfdna_wgs_multiqc_data/multiqc_fastqc.txt",
        qcdir + "/cfdna_wgs_multiqc_data/multiqc_samtools_stats.txt",
        qcdir + "/cfdna_wgs_multiqc_data/multiqc_picard_wgsmetrics.txt",
        qcdir + "/cfdna_wgs_multiqc_data/multiqc_samtools_flagstat.txt",
    params:
        out_dir = qcdir,
        out_name = "cfdna_wgs_multiqc",
        script = cfdna_wgs_scriptdir + "/multiqc.sh",
        threads = cfdna_wgs_threads,
    shell:
        """
        {params.script} \
        "{input}" \
        {params.out_name} \
        {params.out_dir} &> {log}
        """

# Make a tab-separated aggregate QC table
checkpoint cfdna_wgs_make_qc_tsv:
    benchmark: benchdir + "/cfdna_wgs_make_qc_tsv.benchmark.txt",
    input:
        fq = qcdir + "/cfdna_wgs_multiqc_data/multiqc_fastqc.txt",
        mqsam = qcdir + "/cfdna_wgs_multiqc_data/multiqc_samtools_stats.txt",
        mqflag = qcdir + "/cfdna_wgs_multiqc_data/multiqc_samtools_flagstat.txt",
        picard = qcdir + "/cfdna_wgs_multiqc_data/multiqc_picard_wgsmetrics.txt",
        deeptools_frag = qcdir + "/deeptools_frag_lengths.txt",
        deeptools_cov = qcdir + "/cfdna_wgs_coverage.tsv",
    log: logdir + "/cfdna_wgs_make_qc_tsv.log",
    output:
        readqc = qcdir + "/cfdna_wgs_read_qc.tsv",
        fraglen = qcdir + "/cfdna_wgs_frag_len.tsv",
    params:
        script = cfdna_wgs_scriptdir + "/make_qc_tsv.R",
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
