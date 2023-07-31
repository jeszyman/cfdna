#########1#########2#########3#########4#########5#########6#########7#########8
###                                                                          ###
###                    Basic Read Processing of WGS cfDNA                    ###
###                                                                          ###
#########1#########2#########3#########4#########5#########6#########7#########8

rule cfdna_wgs_index:
    output:
        f"{ref_dir}/{{build}}_bwa/{{build}}.fa.sa",
        f"{ref_dir}/{{build}}_bwa/{{build}}.fa",
    params:
        ref_dir = f"{ref_dir}",
        ftp = lambda wildcards: genome_build_map[wildcards.build]['ftp'],
        out_prefix = f"{ref_dir}/{{build}}_bwa/{{build}}.fa",
    shell:
        """
        base=$(basename "{params.ftp}")
        wget -N -P {params.ref_dir} {params.ftp}
        gunzip -c {params.ref_dir}/$base > {params.out_prefix}
        bwa index {params.out_prefix}
        """

rule cfdna_wgs_fastp:
    input:
        read1 = f"{cfdna_wgs_dir}/fastqs/{{library}}_raw_R1.fastq.gz",
        read2 = f"{cfdna_wgs_dir}/fastqs/{{library}}_raw_R2.fastq.gz",
    log:
        html = f"{log_dir}/{{library}}_cfdna_wgs_fastp.html",
    output:
        read1 = f"{cfdna_wgs_dir}/fastqs/{{library}}_proc_R1.fastq.gz",
        read2 = f"{cfdna_wgs_dir}/fastqs/{{library}}_proc_R2.fastq.gz",
        failed = f"{cfdna_wgs_dir}/fastqs/{{library}}_failed_fastp.fastq.gz",
        unpaired1 = f"{cfdna_wgs_dir}/fastqs/{{library}}_unpaired_R1.fastq.gz",
        unpaired2 = f"{cfdna_wgs_dir}/fastqs/{{library}}_unpaired_R2.fastq.gz",
        json = f"{qc_dir}/{{library}}_cfdna_wgs_fastp.json",
        cmd = f"{qc_dir}/{{library}}_cfdna_wgs_fastp.log",
    params:
        script = f"{cfdna_script_dir}/fastp.sh",
        threads = threads,
    resources:
        mem_mb = 500,
    shell:
        """
        {params.script} \
        {input.read1} \
        {input.read2} \
        {log.html} \
        {output.json} \
        {output.read1} \
        {output.read2} \
        {output.failed} \
        {output.unpaired1} \
        {output.unpaired2} \
        {params.threads} &> {output.cmd}
        """

# Align reads with BWA
rule frag_align:
    benchmark: f"{bench_dir}/{{library}}_{{build}}_frag_align.benchmark.txt",
    input:
        ref = f"{ref_dir}/{{build}}_bwa/{{build}}.fa.sa",
        read1 = f"{cfdna_wgs_dir}/fastqs/{{library}}_proc_R1.fastq.gz",
        read2 = f"{cfdna_wgs_dir}/fastqs/{{library}}_proc_R2.fastq.gz",
    log: f"{log_dir}/{{library}}_{{build}}frag_align.log",
    output:
        sort = f"{cfdna_wgs_dir}/bams/{{library}}_{{build}}_raw.bam",
        index = f"{cfdna_wgs_dir}/bams/{{library}}_{{build}}_raw.bam.bai",
    params:
        ref = f"{ref_dir}/{{build}}_bwa/{{build}}.fa",
        script = f"{cfdna_script_dir}/align.sh",
        threads = 4,
    resources:
        mem_mb = 500,
    shell:
        """
        {params.script} \
        {params.ref} \
        {input.read1} \
        {input.read2} \
        {params.threads} \
        {output.sort} &> {log}
        """

# Remove PCR duplicates from aligned reads
rule frag_dedup:
    benchmark: f"{bench_dir}/{{library}}_{{build}}_frag_dedup.benchmark.txt",
    input: f"{cfdna_wgs_dir}/bams/{{library}}_{{build}}_raw.bam",
    log: f"{log_dir}/{{library}}_{{build}}_frag_dedup.log",
    output: f"{cfdna_wgs_dir}/bams/{{library}}_{{build}}_dedup.bam",
    params:
        script = f"{cfdna_script_dir}/dedup.sh",
        threads = threads,
    shell:
        """
        {params.script} \
        {input} \
        {output} \
        {params.threads} &> {log}
        """

rule frag_filter_alignment:
    benchmark: f"{bench_dir}/{{library}}_{{build}}_frag_filter_alignment.benchmark.txt",
    input: f"{cfdna_wgs_dir}/bams/{{library}}_{{build}}_dedup.bam",
    log: f"{log_dir}/{{library}}_{{build}}_frag_filter_alignment.log",
    output: f"{cfdna_wgs_dir}/bams/{{library}}_{{build}}_filt.bam",
    params:
        script = f"{cfdna_script_dir}/filter_alignment.sh",
        threads = threads,
    shell:
        """
        {params.script} \
        {input} \
        {params.threads} \
        {output} &> {log}
        """

# Get read quality by FASTQC
rule frag_fastqc:
    benchmark: f"{bench_dir}/{{library}}_{{processing}}_{{read}}_frag_fastqc.benchmark.txt",
    input: f"{cfdna_wgs_dir}/fastqs/{{library}}_{{processing}}_{{read}}.fastq.gz",
    log: f"{log_dir}/{{library}}_{{processing}}_{{read}}_frag_fastqc.log",
    output:
        f"{qc_dir}/{{library}}_{{processing}}_{{read}}_fastqc.html",
        f"{qc_dir}/{{library}}_{{processing}}_{{read}}_fastqc.zip",
    params:
        outdir = f"{qc_dir}",
        script = f"{cfdna_script_dir}/fastqc.sh",
        threads = threads,
    shell:
        """
        {params.script} \
        {input} \
        {params.outdir} \
        {params.threads} &> {log}
        """

# Get alignment QC using samtools
rule frag_alignment_qc:
    input: f"{cfdna_wgs_dir}/bams/{{library}}_{{build}}_{{processing}}.bam",
    log:
        flagstat = f"{log_dir}/{{library}}_{{build}}_{{processing}}_flagstat_frag_alignment_qc.log",
        samstat = f"{log_dir}/{{library}}_{{build}}_{{processing}}_samstats_frag_alignment_qc.log",
    output:
        flagstat = f"{qc_dir}/{{library}}_{{build}}_{{processing}}_flagstat.txt",
        samstat = f"{qc_dir}/{{library}}_{{build}}_{{processing}}_samstats.txt",
    params:
        script = f"{cfdna_script_dir}/alignment_qc.sh",
        threads = threads,
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
    benchmark: f"{bench_dir}/{{library}}_{{build}}_frag_picard_depth.benchmark.txt",
    input:
        ref = f"{ref_dir}/{{build}}_bwa/{{build}}.fa",
        bam = f"{cfdna_wgs_dir}/bams/{{library}}_{{build}}_filt.bam",
    log: f"{cfdna_wgs_dir}/bams/{{library}}_{{build}}_frag_picard_depth.log",
    output: f"{qc_dir}/{{library}}_{{build}}_picard_depth.txt",
    params:
        script = f"{cfdna_script_dir}/picard_depth.sh",
        threads = threads,
    shell:
        """
        {params.script} \
        {input.bam} \
        {input.ref} \
        {output}
        """

# Get fragment sizes using deepTools
rule frag_bampefragsize:
    conda: "deeptools",
    input:
        lambda wildcards: expand(f"{cfdna_wgs_dir}/bams/{{library}}_{{build}}_filt.bam",
                                 library = lib_map[wildcards.lib_set]['libs'],
                                 build = lib_map[wildcards.lib_set]['build']),
    log: f"{log_dir}/{{lib_set}}_bampefragsize.log",
    output:
        raw = f"{qc_dir}/deeptools_{{lib_set}}_lengths.txt",
        hist = f"{qc_dir}/deeptools_{{lib_set}}_lengths.png",
    params:
        blacklist = lambda wildcards: lib_map[wildcards.lib_set]['blacklist'],
        script = f"{cfdna_script_dir}/bampefragsize.sh",
        threads = threads,
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

# Make deepTools plotCoverage coverage maps for all filtered bams
rule frag_plotcoverage:
    conda: "deeptools",
    input:
        lambda wildcards: expand(f"{cfdna_wgs_dir}/bams/{{library}}_{{build}}_filt.bam",
                                 library = lib_map[wildcards.lib_set]['libs'],
                                 build = lib_map[wildcards.lib_set]['build']),
    log: f"{log_dir}/{{lib_set}}_frag_plotcoverage.log",
    output:
        raw = f"{qc_dir}/{{lib_set}}_frag_coverage.tsv",
        plot = f"{qc_dir}/{{lib_set}}_frag_coverage.pdf",
    params:
        blacklist = lambda wildcards: lib_map[wildcards.lib_set]['blacklist'],
        script = f"{cfdna_script_dir}/plotcoverage.sh",
        threads = threads,
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
    input:
        lambda wildcards: expand(f"{qc_dir}/{{library}}_cfdna_wgs_fastp.json",
                                 library = lib_map[wildcards.lib_set]['libs']),
        lambda wildcards: expand(f"{qc_dir}/{{library}}_{{processing}}_{{read}}_fastqc.zip",
                                 library = lib_map[wildcards.lib_set]['libs'],
                                 processing = ['raw','proc',],
                                 read = ['R1','R2']),
        lambda wildcards: expand(f"{qc_dir}/{{library}}_{{build}}_{{processing}}_samstats.txt",
                                 library = lib_map[wildcards.lib_set]['libs'],
                                 build = lib_map[wildcards.lib_set]['build'],
                                 processing = ['raw','dedup','filt']),
        lambda wildcards: expand(f"{qc_dir}/{{library}}_{{build}}_{{processing}}_flagstat.txt",
                                 library = lib_map[wildcards.lib_set]['libs'],
                                 build = lib_map[wildcards.lib_set]['build'],
                                 processing = ['raw','dedup','filt']),
        lambda wildcards: expand(f"{qc_dir}/{{library}}_{{build}}_picard_depth.txt",
                                 library = lib_map[wildcards.lib_set]['libs'],
                                 build = lib_map[wildcards.lib_set]['build']),
        f"{qc_dir}/deeptools_{{lib_set}}_lengths.txt",
        f"{qc_dir}/{{lib_set}}_frag_coverage.tsv",
    log: f"{log_dir}/{{lib_set}}_cfdna_wgs_multiqc.log",
    output:
        f"{qc_dir}/{{lib_set}}_cfdna_wgs_multiqc.html",
        f"{qc_dir}/{{lib_set}}_cfdna_wgs_multiqc_fastqc.txt",
        f"{qc_dir}/{{lib_set}}_cfdna_wgs_multiqc_data/multiqc_samtools_stats.txt",
        f"{qc_dir}/{{lib_set}}_cfdna_wgs_multiqc_data/multiqc_picard_wgsmetrics.txt",
        f"{qc_dir}/{{lib_set}}_cfdna_wgs_multiqc_data/multiqc_samtools_flagstat.txt",
    params:
        out_dir = f"{qc_dir}",
        out_name = "frag_multiqc",
        script = f"{cfdna_script_dir}/multiqc.sh",
        threads = threads,
    shell:
        """
        {params.script} \
        "{input}" \
        {params.out_name} \
        {params.out_dir} &> {log}
        """

# Make a tab-separated aggregate QC table
rule make_cfdna_wgs_qc_tsv:
    input:
        fq = f"{qc_dir}/{{lib_set}}_cfdna_wgs_multiqc_data/multiqc_fastqc.txt",
        mqsam = f"{qc_dir}/{{lib_set}}_cfdna_wgs_multiqc_data/multiqc_samtools_stats.txt",
        mqflag = f"{qc_dir}/{{lib_set}}_cfdna_multiqc_data/multiqc_samtools_flagstat.txt",
        picard = f"{qc_dir}/{{lib_set}}_multiqc_data/multiqc_picard_wgsmetrics.txt",
        deeptools_frag = f"{qc_dir}/{{lib_set}}_deeptools_cfdna_wgs_lengths.txt",
        deeptools_cov = f"{qc_dir}/{{lib_set}}_cfdna_wgs_coverage.tsv",
    log: f"{log_dir}/{{lib_set}}_cfdna_wgs_make_qc_tsv.log",
    output:
        readqc = f"{qc_dir}/{{lib_set}}_cfdna_wgs_read_qc.tsv",
        fraglen = f"{qc_dir}/{{lib_set}}_cfdna_wgs_len.tsv",
    params:
        script = f"{cfdna_script_dir}/make_qc_tsv.R",
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
