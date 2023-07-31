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
