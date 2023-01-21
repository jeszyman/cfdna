#########1#########2#########3#########4#########5#########6#########7#########8
#                                                                              #
#                    Basic Read and Alignment Processing of                    #
#                    Cell-free DNA Whole Genome Sequencing                     #
#                                                                              #
#########1#########2#########3#########4#########5#########6#########7#########8

# Make alignment index
#  Note: Upon first run, this rule will touch an empty file with the same path
#        as the input fasta. Thereafter, you can avoid repeat indexing when the
#        rule "sees" this empty file. For repo intergration testing with an
#        external reference, indexing can likewise be avoided with this empty
#        file at the external index location.

rule cfdna_wgs_index:
    benchmark: benchdir + "/cfdna_wgs_index.benchmark.txt",
    container: cfdna_wgs_container,
    input: genome_fasta,
    log: logdir + "/cfdna_wgs_index.log",
    output: done = touch(genome_ref),
    params:
        out_prefix = genome_ref,
        script = cfdna_wgs_scriptdir + "/index.sh",
        threads = cfdna_wgs_threads,
    shell:
        """
        bwa index -p {params.out_prefix} {input} &> {log}
        """

# Make a file of blacklist-filtered autosomal regions
rule cfdna_wgs_make_keep_bed:
    benchmark: benchdir + "/cfdna_wgs_make_keep_bed.benchmark.txt",
    container: cfdna_wgs_container,
    input:
        blacklist = blklist,
        chrom_sizes = chrom_sizes,
    log: logdir + "/cfdna_wgs_make_keep_bed.log",
    output:
        autosome_bed = autosome_bed,
        keep_bed = keep_bed,
    params:
        script = cfdna_wgs_scriptdir + "/make_keep_bed.sh",
    shell:
        """
        {params.script} \
        {input.blacklist} \
        {input.chrom_sizes} \
        {output.autosome_bed} \
        {output.keep_bed} &> {log}
        """

# Adapter-trim and QC reads with fastp
rule cfdna_wgs_fastp:
    benchmark: benchdir + "/{library}_cfdna_wgs_fastp.benchmark.txt",
    container: cfdna_wgs_container,
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
    container: cfdna_wgs_container,
    input:
        ref = genome_ref,
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
    container: cfdna_wgs_container,
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
    container: cfdna_wgs_container,
    input:
        bam = cfdna_wgs_bams + "/{library}_dedup.bam",
        keep_bed = keep_bed,
    log: logdir + "/{library}_cfdna_wgs_filter_alignment.log",
    output: cfdna_wgs_bams + "/{library}_filt.bam",
    params:
        script = cfdna_wgs_scriptdir + "/filter_alignment.sh",
        threads = cfdna_wgs_threads,
    shell:
        """
        {params.script} \
        {input.bam} \
        {input.keep_bed} \
        {params.threads} \
        {output} &> {log}
        """

# Get read quality by FASTQC
rule cfdna_wgs_fastqc:
    benchmark: benchdir+ "/{library}_{processing}_{read}_cfdna_wgs_fastqc.benchmark.txt",
    container: cfdna_wgs_container,
    input: cfdna_wgs_fastqs + "/{library}_{processing}_{read}.fastq.gz",
    log: logdir + "/{library}_{processing}_{read}_cfdna_wgs_fastqc.log",
    output:
        qc + "/{library}_{processing}_{read}_fastqc.html",
        qc + "/{library}_{processing}_{read}_fastqc.zip",
    params:
        outdir = qc,
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
    container: cfdna_wgs_container,
    input: cfdna_wgs_bams + "/{library}_{processing}.bam",
    log:
        flagstat = logdir + "/{library}_{processing}_flagstat_cfdna_wgs_alignment_qc.log",
        samstat = logdir + "/{library}_{processing}_samstats_cfdna_wgs_alignment_qc.log",
    output:
        flagstat = qc + "/{library}_{processing}_flagstat.txt",
        samstat = qc + "/{library}_{processing}_samstats.txt",
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
    container: cfdna_wgs_container,
    input: cfdna_wgs_bams + "/{library}_filt.bam",
    log: logdir + "/{library}_cfdna_wgs_picard_depth.log",
    output: qc + "/{library}_picard_depth.txt",
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
    container: cfdna_wgs_container,
    input: expand(cfdna_wgs_bams + "/{library}_filt.bam", library = CFDNA_WGS_LIBRARIES),
    log: logdir + "/cfdna_wgs_bampefragsize.log",
    output:
        raw = qc + "/deeptools_frag_lengths.txt",
        hist = qc + "/deeptools_frag_lengths.png",
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
    container: cfdna_wgs_container,
    input: cfdna_wgs_bams + "/{library}_filt.bam",
    log: logdir + "/{library}_cfdna_wgs_bamcoverage.log",
    output: qc + "/{library}_bamcoverage.bg",
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
    container: cfdna_wgs_container,
    input: expand(cfdna_wgs_bams + "/{library}_filt.bam", library = CFDNA_WGS_LIBRARIES),
    log: logdir + "/cfdna_wgs_plotcoverage.log",
    output:
        raw = qc + "/cfdna_wgs_coverage.tsv",
        plot = qc + "/cfdna_wgs_coverage.pdf",
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
    container: cfdna_wgs_container,
    input:
        expand(logdir + "/{library}_cfdna_wgs_fastp.json", library = CFDNA_WGS_LIBRARIES),
        expand(qc + "/{library}_{processing}_{read}_fastqc.zip", library = CFDNA_WGS_LIBRARIES, processing = ["raw", "processed", "unpaired"], read = ["R1","R2"]),
        expand(qc + "/{library}_{processing}_samstats.txt", library = CFDNA_WGS_LIBRARIES, processing = ["raw","filt"]),
        expand(qc + "/{library}_{processing}_flagstat.txt", library = CFDNA_WGS_LIBRARIES, processing = ["raw","filt"]),
        expand(qc + "/{library}_picard_depth.txt", library = CFDNA_WGS_LIBRARIES),
        qc + "/deeptools_frag_lengths.txt",
        qc + "/cfdna_wgs_coverage.tsv",
    log: logdir + "/cfdna_wgs_multiqc.log",
    output:
        qc + "/cfdna_wgs_multiqc.html",
        qc + "/cfdna_wgs_multiqc_data/multiqc_fastqc.txt",
        qc + "/cfdna_wgs_multiqc_data/multiqc_samtools_stats.txt",
        qc + "/cfdna_wgs_multiqc_data/multiqc_picard_wgsmetrics.txt",
        qc + "/cfdna_wgs_multiqc_data/multiqc_samtools_flagstat.txt",
    params:
        out_dir = qc,
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
    container: cfdna_wgs_container,
    input:
        fq = qc + "/cfdna_wgs_multiqc_data/multiqc_fastqc.txt",
        mqsam = qc + "/cfdna_wgs_multiqc_data/multiqc_samtools_stats.txt",
        mqflag = qc + "/cfdna_wgs_multiqc_data/multiqc_samtools_flagstat.txt",
        picard = qc + "/cfdna_wgs_multiqc_data/multiqc_picard_wgsmetrics.txt",
        deeptools_frag = qc + "/deeptools_frag_lengths.txt",
        deeptools_cov = qc + "/cfdna_wgs_coverage.tsv",
    log: logdir + "/cfdna_wgs_make_qc_tsv.log",
    output:
        readqc = qc + "/cfdna_wgs_read_qc.tsv",
        fraglen = qc + "/cfdna_wgs_frag_len.tsv",
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
