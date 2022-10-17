##############################
###   cfDNA WGS Pipeline   ###
##############################

# Raw read processing with fastp
rule fastp:
    container:
        config["container"]["default"],
    input:
        read1 = cfdna_wgs_fastqs + "/{library}_raw_R1.fastq.gz",
        read2 = cfdna_wgs_fastqs + "/{library}_raw_R2.fastq.gz",
    log:
        cmd = config["logdir"] + "/{library}_fastp.log",
        html = config["logdir"] + "/{library}_fastp.html",
        json = config["logdir"] + "/{library}_fastp.json",
    output:
        read1 = cfdna_wgs_fastqs + "/{library}_processed_R1.fastq.gz",
        read2 = cfdna_wgs_fastqs + "/{library}_processed_R2.fastq.gz",
        failed = cfdna_wgs_fastqs + "/{library}_failed_fastp.fastq.gz",
        unpaired1 = cfdna_wgs_fastqs + "/{library}_unpaired_R1.fastq.gz",
        unpaired2 = cfdna_wgs_fastqs + "/{library}_unpaired_R2.fastq.gz",
    params:
        script = config["scriptdir"]["cfdna_wgs"] + "/fastp.sh",
        threads = config["threads"]["default"],
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

rule index:
    container:
        config["container"]["cfdna_wgs"],
    input:
        config["genome_fasta"],
    output:
        done = touch(genome_ref)
    params:
        out_prefix = genome_ref
    shell:
        """
        bwa index -p {params.out_prefix} {input}
        """

# BWA alignment
rule align:
    benchmark:
        config["logdir"] + "/{library}_align.benchmark.txt",
    container:
        config["container"]["default"],
    input:
        ref = genome_ref,
        r1 = cfdna_wgs_fastqs + "/{library}_processed_R1.fastq.gz",
        r2 = cfdna_wgs_fastqs + "/{library}_processed_R2.fastq.gz",
    log:
        config["logdir"] + "/{library}_align.log",
    output:
        sort = cfdna_wgs_bams + "/{library}_raw.bam",
        index = cfdna_wgs_bams + "/{library}_raw.bam.bai",
    params:
        script = config["scriptdir"]["cfdna_wgs"] + "/align.sh",
        threads = config["threads"]["bwa"]
    shell:
        """
        {params.script} \
        {input.ref} \
        {input.r1} \
        {input.r2} \
        {params.threads} \
        {output.sort} &> {log}
	"""

# Remove PCR duplicates from aligned reads
rule dedup:
    container:
        config["container"]["cfdna_wgs"],
    input:
        cfdna_wgs_bams + "/{library}_raw.bam",
    log:
        config["logdir"] + "/{library}_cfdna_wgs_bam_dedup.log",
    output:
        cfdna_wgs_bams + "/{library}_dedup.bam",
    params:
        script = config["scriptdir"]["cfdna_wgs"] + "/dedup.sh",
        threads = config["threads"]["bwa"],
    shell:
        """
        {params.script} \
        {input} \
        {output} \
        {threads} &> {log}
        """

# Removes unmapped, not primary, and duplicate reads. Additionally, quality filters by config variable.
rule alignment_filtering:
    container:
        config["container"]["cfdna_wgs"],
    input:
        cfdna_wgs_bams + "/{library}_dedup.bam",
    log:
        config["logdir"] + "/{library}_alignment_filtering.log",
    output:
        cfdna_wgs_bams + "/{library}_filt.bam",
    params:
        keepbed = config["keepbed"],
        script = config["scriptdir"]["cfdna_wgs"] + "/alignment_filtering.sh",
        threads = config["threads"]["default"],
    shell:
        """
        {params.script} \
        {input} \
        {params.keepbed} \
        {params.threads} \
        {output} &> {log}
        """

# FastQC
rule fastqc:
    container:
        config["container"]["default"]
    input:
        cfdna_wgs_fastqs + "/{library}_{processing}_{read}.fastq.gz",
    log:
        config["logdir"] + "/{library}_{processing}_{read}_fastqc.log",
    output:
        config["qcdir"] + "/{library}_{processing}_{read}_fastqc.html",
        config["qcdir"] + "/{library}_{processing}_{read}_fastqc.zip",
    params:
        outdir = config["qcdir"],
        script = config["scriptdir"]["cfdna_wgs"] + "/fastqc_wrapper.sh",
        threads = config["threads"]["default"],
    shell:
        """
        {params.script} \
        {input} \
        {params.outdir} \
        {params.threads} &> {log}
        """

# Alignment samtools QC
rule alignment_qc:
    container:
        config["container"]["cfdna_wgs"],
    input:
        cfdna_wgs_bams + "/{library}_{processing}.bam",
    log:
        flagstat = config["logdir"] + "/{library}_{processing}_flagstat.log",
        samstat = config["logdir"] + "/{library}_{processing}_samstat.log",
    output:
        flagstat = config["qcdir"] + "/{library}_{processing}_flagstat.txt",
        samstat = config["qcdir"] + "/{library}_{processing}_samstats.txt",
    params:
        script = config["scriptdir"]["cfdna_wgs"] + "/alignment_qc.sh",
        threads = config["threads"]["default"],
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

# Sequencing depth via Picard
rule picard_wgs:
    container:
        config["container"]["cfdna_wgs"],
    input:
        cfdna_wgs_bams + "/{library}_filt.bam",
    log:
        config["logdir"] + "/{library}_picard_wgs.log",
    output:
        config["qcdir"] + "/{library}_picard_wgs.txt",
    params:
        script = config["scriptdir"]["cfdna_wgs"] + "/CollectWgsMetrics_wrapper.sh",
    shell:
        """
        {params.script} \
        {input} \
        {config[picard_jar]} \
        {config[genome_fasta]} \
        {output}
        """

# Fragment sizes by deepTools
rule deeptools_bampefragmentsize:
    container:
        config["container"]["cfdna_wgs"],
    input:
        expand(cfdna_wgs_bams + "/{library}_filt.bam", library = CFDNA_WGS_LIBRARIES),
    log:
        config["logdir"] + "/bampefragmentsize.txt",
    output:
        hist = config["qcdir"] + "/deeptools_frag_lengths.png",
        raw = config["qcdir"] + "/deeptools_frag_lengths.txt",
    params:
        blacklist = config["blacklist"],
        script = config["scriptdir"]["cfdna_wgs"] + "/bamPEFragmentSize_wrapper.sh",
        threads = config["threads"]["default"],
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
rule bamcoverage:
    container:
        config["container"]["cfdna_wgs"],
    input:
        cfdna_wgs_bams + "/{library}_filt.bam",
    log:
        config["logdir"] + "/{library}_bamcoverage.log",
    output:
        config["qcdir"] + "/{library}_bamcoverage.bg",
    params:
        bin = "10000",
        blacklist = config["blacklist"],
        script = config["scriptdir"]["cfdna_wgs"] + "/bamcoverage.sh",
        threads = config["threads"]["default"],
    shell:
        """
        {params.script} \
        {input} \
        {output} \
        {params.bin} \
        {params.blacklist} \
        {params.threads} &> {log}
        """

# deeptools plotCoverage on all filtered bams
rule plot_coverage:
    container:
        config["container"]["cfdna_wgs"],
    input:
        expand(cfdna_wgs_bams + "/{library}_filt.bam", library = CFDNA_WGS_LIBRARIES),
    log:
        config["logdir"] + "/plot_coverage.log",
    output:
        raw = config["qcdir"] + "/coverage.tsv",
        plot = config["qcdir"] + "/coverage.pdf",
    params:
        blacklist = config["blacklist"],
        script = config["scriptdir"]["cfdna_wgs"] + "/plot_coverage.sh",
    shell:
        """
        {params.script} \
        "{input}" \
        {params.blacklist} \
        {config[threads][default]} \
        {output.raw} \
        {output.plot} &> {log}
        """

rule cfdna_wgs_multiqc:
    container:
        config["container"]["cfdna_wgs"],
    input:
        expand(config["logdir"] + "/{library}_fastp.json",
            library = CFDNA_WGS_LIBRARIES),
        expand(config["qcdir"] + "/{library}_{processing}_{read}_fastqc.zip",
            library = CFDNA_WGS_LIBRARIES,
            processing = ["raw", "processed", "unpaired"],
            read = ["R1","R2"]),
        expand(config["qcdir"] + "/{library}_{processing}_flagstat.txt",
            library = CFDNA_WGS_LIBRARIES,
            processing = ["raw", "dedup", "filt"]),
        expand(config["qcdir"] + "/{library}_{processing}_samstats.txt",
            library = CFDNA_WGS_LIBRARIES,
            processing = ["raw", "dedup", "filt"]),
        expand(config["qcdir"] + "/{library}_picard_wgs.txt",
            library = CFDNA_WGS_LIBRARIES),
        config["qcdir"] + "/deeptools_frag_lengths.txt",
        config["qcdir"] + "/coverage.tsv",
    log:
        config["logdir"] + "/cfdna_wgs_multiqc.log"
    output:
        config["qcdir"] + "/all_cfdna_wgs.html",
    params:
        out_dir = config["qcdir"],
        out_name = "all_cfdna_wgs",
        script = config["scriptdir"]["cfdna_wgs"] + "/cfdna_wgs_multiqc.sh",
    shell:
        """
        {params.script} \
        "{input}" \
        {params.out_name} \
        {params.out_dir} &> {log}
        """

#  Notes:
#  This makes an aggregate table of QC values. The subsequent downsampling
#  step only runs if read numbers are above a certain threshold. See also
#  the int_test.smk for function using this output table.

checkpoint make_qc_tbl:
    container:
        config["container"]["cfdna_wgs"],
    input:
        fq = config["qcdir"] + "/all_cfdna_wgs_data/multiqc_fastqc.txt",
        sam = config["qcdir"] + "/all_cfdna_wgs_data/multiqc_samtools_stats.txt",
        flag = config["qcdir"] + "/all_cfdna_wgs_data/multiqc_samtools_flagstat.txt",
	picard = config["qcdir"] + "/all_cfdna_wgs_data/multiqc_picard_wgsmetrics.txt",
        deeptools_frag = config["qcdir"] + "/deeptools_frag_lengths.txt",
        deeptools_cov = config["qcdir"] + "/coverage.tsv"
    log:
        config["logdir"] + "/read_qc.log"
    output:
        readqc = config["qcdir"] + "/cfdna_wgs_read_qc.tsv",
        fraglen = config["qcdir"] + "/cfdna_wgs_frag_len.tsv",
    params:
        script = config["scriptdir"]["cfdna_wgs"] + "/make_qc_tbl.R"
    shell:
        """
        Rscript {params.script} \
        {input.fq} \
        {input.sam} \
        {input.flag} \
        {input.picard} \
        {input.deeptools_frag} \
        {input.deeptools_cov} \
        {output.readqc} \
        {output.fraglen} >& {log}
        """
