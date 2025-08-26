###############################
###   cfDNA WGS Snakefile   ###
###############################

#########1#########2#########3#########4#########5#########6#########7#########8
# A snakefile for basic processing of cell-free DNA whole-genome sequecing data.

# ---   Dependencies   --- #
# ------------------------ #

# ./config/cfdna-wgs-conda-env.yaml, a conda environment file
# Scripts within ./scripts

# ---   Configuration Parameters   --- #
# ------------------------------------ #

# Parameters to be defined at the configuration yaml include:
#
# available_concurrency:
# data-temp-dir:
# mosdepth-quant-levels:
#
# ref_assemblies:
#   <ASSEMBLY ID>:
#     url:
#     name:
#     input:

# Example minimal configuration yaml:
# This is a nested map, e.g.:
# ref_assemblies:
#   ncbi_decoy_hg38:
#     url: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
#     name: ncbi_decoy_hg38
#     input: GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fa.gz


# {Parameters} to be defined at the level of the snakefile include
#
# data_dir
# config['cfdna-scripts-dir']

# {{Wildcards}} to be defined in the wrapper snakefile rule all include:
#
# library_id
# processing
# read
# ref_name
# align_method

#########1#########2#########3#########4#########5#########6#########7#########8
#
# This is a modular snakefile, intended to be incorporated into a larger
# workflow using the "include:" directive. (See
# https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html)
#
#########1#########2#########3#########4#########5#########6#########7#########8




#########1#########2#########3#########4#########5#########6#########7#########8 
#
# This snakefile uses a common data directory structure:
#
# |-- <MAIN PROJECT DATA DIRECTORY>, e.g. lung_cancer, variable {data_dir}>
#     |-- inputs
#     |-- ref
#     |-- logs
#     |-- <ANALYSIS LABEL>, e.g. emseq, example subdirectories: 
#         |-- qc
#         |-- fastqs
#         |-- bams
#

# Specifically here, most ouputs will return to {data_dir}/cfdna-wgs

#nil

#nil
rule cfdna_wgs_fastp:
    #
    # fastp for cfDNA WGS. Uses a set thread count of 8. Adapters are
    # auto-detected.
    #
    conda:
        "../config/cfdna-wgs-conda-env.yaml"
    input:
        r1 = f"{config['cfdna-wgs-dir']}/fastqs/{{library_id}}.raw_R1.fastq.gz",
        r2 = f"{config['cfdna-wgs-dir']}/fastqs/{{library_id}}.raw_R2.fastq.gz",
    log:
        html = f"{config['log-dir']}/{{library_id}}_cfdna_wgs_fastp.html",
        json = f"{config['log-dir']}/{{library_id}}_cfdna_wgs_fastp.json",
        run = f"{config['log-dir']}/{{library_id}}_cfdna_wgs_fastp.log",
    output:
        failed = f"{config['cfdna-wgs-dir']}/fastqs/{{library_id}}.failed.fastq.gz",
        r1 = f"{config['cfdna-wgs-dir']}/fastqs/{{library_id}}.processed_R1.fastq.gz",
        r2 = f"{config['cfdna-wgs-dir']}/fastqs/{{library_id}}.processed_R2.fastq.gz",
        up1 = f"{config['cfdna-wgs-dir']}/fastqs/{{library_id}}.unpaired_R1.fastq.gz",
        up2 = f"{config['cfdna-wgs-dir']}/fastqs/{{library_id}}.unpaired_R2.fastq.gz",
    resources:
        concurrency=12,
    threads:
        10
    shell:
        """
        fastp --detect_adapter_for_pe \
        --failed_out {output.failed} \
        --in1 {input.r1} --in2 {input.r2} \
        --html {log.html} --json {log.json} \
        --out1 {output.r1} --out2 {output.r2} \
        --unpaired1 {output.up1} --unpaired2 {output.up2} \
        --thread {threads} &> {log.run}
        """
rule cfdna_wgs_fastqc:
    conda:
        "../config/cfdna-wgs-conda-env.yaml"
    input:
        f"{config['cfdna-wgs-dir']}/fastqs/{{library_id}}.{{processing}}_{{read}}.fastq.gz",
    log:
        f"{config['log-dir']}/{{library_id}}.{{processing}}_{{read}}_cfdna_wgs_fastqc.log",
    output:
        f"{config['cfdna-wgs-dir']}/qc/{{library_id}}.{{processing}}_{{read}}_fastqc.html",
        f"{config['cfdna-wgs-dir']}/qc/{{library_id}}.{{processing}}_{{read}}_fastqc.zip",
    params:
        outdir = f"{config['cfdna-wgs-dir']}/qc",
    resources:
        concurrency = 25,
    threads:
        2
    shell:
        """
        fastqc \
        --outdir {params.outdir} \
        --quiet \
        --svg \
        --threads {threads} \
        {input} &> {log}
        """
#########1#########2#########3#########4#########5#########6#########7#########8
rule cfdna_wgs_bwa_index:
    #
    # Builds BWA reference off of an existing fasta file in the
    # {data_dir}/inputs directory. Uses a nested map from the config yaml. Also
    # created are a fasta index (.fai) and a bed file of the primary assembly
    # (numbered chromosomes and X and Y, no mitochondria or other contigs).
    #
    conda:
        "../config/cfdna-wgs-conda-env.yaml"
    input:
        lambda wildcards: f"{config['data-dir']}/inputs/{config['ref_assemblies'][wildcards.name]['input']}",
    output:
        fa = f"{config['data-dir']}/ref/bwa/{{name}}/{{name}}.fa",
        fai = f"{config['data-dir']}/ref/bwa/{{name}}/{{name}}.fa.fai",
        bed = f"{config['data-dir']}/ref/bwa/{{name}}/{{name}}.primary.bed",
        amb = f"{config['data-dir']}/ref/bwa/{{name}}/{{name}}.amb",
        ann = f"{config['data-dir']}/ref/bwa/{{name}}/{{name}}.ann",
        bwt = f"{config['data-dir']}/ref/bwa/{{name}}/{{name}}.bwt",
        pac = f"{config['data-dir']}/ref/bwa/{{name}}/{{name}}.pac",
        sa  = f"{config['data-dir']}/ref/bwa/{{name}}/{{name}}.sa",
    params:
        bwa_prefix = lambda wildcards: f"{config['data-dir']}/ref/bwa/{wildcards.name}/{wildcards.name}",
        script = f"{config['cfdna-scripts-dir']}/cfdna_wgs_bwa_index.sh",
    log:
        f"{config['log-dir']}/{{name}}_bwa_index.log"
    shell:
        """
        {params.script} \
        {input} \
        {output.fa} \
        {output.fai} \
        {output.bed} \
        {params.bwa_prefix} \
        {log}
        """
rule cfdna_wgs_bwa_mem:
    conda:
        "../config/cfdna-wgs-conda-env.yaml"
    input:
        r1 = f"{config['cfdna-wgs-dir']}/fastqs/{{library_id}}.processed_R1.fastq.gz",
        r2 = f"{config['cfdna-wgs-dir']}/fastqs/{{library_id}}.processed_R2.fastq.gz",
        sa_check = f"{config['data-dir']}/ref/bwa/{{ref_name}}/{{ref_name}}.sa",
    output:
        bam = f"{config['cfdna-wgs-dir']}/bams/{{library_id}}.{{ref_name}}.bwa.coorsort.bam",
    params:
        ref = f"{config['data-dir']}/ref/bwa/{{ref_name}}/{{ref_name}}",
    resources:
        concurrency = 100,
    threads:
        50
    shell:
        """
        bwa mem -M -t {threads} \
        {params.ref} {input.r1} {input.r2} \
        | samtools view -@ 4 -Sb - -o - \
        | samtools sort -@ 4 - -o {output.bam}
        samtools index -@ 4 {output.bam}
        """
rule cfdna_wgs_bam_dedup:
    #
    # 1) Name sort, required by fixmate
    # 2) Fixmate adds mate-pair info needed for deduping
    # 3) Coordinate sorting, required by markdup
    # 4) Markdup REMOVING PCR duplicates
    #
    conda:
        "../config/cfdna-wgs-conda-env.yaml",
    input:
        f"{config['cfdna-wgs-dir']}/bams/{{library_id}}.{{ref_name}}.{{align_method}}.coorsort.bam",
    log:
        f"{config['log-dir']}/{{library_id}}.{{ref_name}}.{{align_method}}_cfdna_wgs_bam_dedup.log",
    output:
        bam = f"{config['cfdna-wgs-dir']}/bams/{{library_id}}.{{ref_name}}.{{align_method}}.dedup.coorsort.bam",
        bai = f"{config['cfdna-wgs-dir']}/bams/{{library_id}}.{{ref_name}}.{{align_method}}.dedup.coorsort.bam.bai",
    params:
        tmp_dir = config["local-tmp-dir"]
    shell:
        """
        rm -rf  {params.tmp_dir}/{wildcards.library_id}.namesort
        mkdir -p {params.tmp_dir}/{wildcards.library_id}.namesort
        samtools sort -@ 8 -n -T {params.tmp_dir}/{wildcards.library_id}.namesort -o - {input} \
        | samtools fixmate -@ 8 -m - - \
        | samtools sort -@ 8 -T {params.tmp_dir}/{wildcards.library_id}.namesort -o - - \
        | samtools markdup -@ 8 -r -T {params.tmp_dir}/{wildcards.library_id}.namesort - {output.bam}
        samtools index -@ 4 {output.bam}
        """
rule cfdna_wgs_bam_filt:
    #
    # Excludes any unmapped (0x4),
    #  not primary alignment (0x100),
    #  or duplicates (0x400)
    #
    # Only MAPQ > 20
    #
    # Restrict to primary chromosomes
    #
    # DO NOT restrict to "proper pairs"- this clips long cfDNA fragments!
    #
    conda:
        "../config/cfdna-wgs-conda-env.yaml",
    input:
        bam = f"{config['cfdna-wgs-dir']}/bams/{{library_id}}.{{ref_name}}.{{align_method}}.dedup.coorsort.bam",
        bed = f"{config['data-dir']}/ref/bwa/{{ref_name}}/{{ref_name}}.primary.bed",
    log:
        f"{config['log-dir']}/{{library_id}}.{{ref_name}}.{{align_method}}_cfdna_wgs_bam_filt.log",
    output:
        bam = f"{config['cfdna-wgs-dir']}/bams/{{library_id}}.{{ref_name}}.{{align_method}}.dedup.coorsort.filt.bam",
        bai = f"{config['cfdna-wgs-dir']}/bams/{{library_id}}.{{ref_name}}.{{align_method}}.dedup.coorsort.filt.bam.bai",
    shell:
        """
        samtools view -@ 8 -b -F 1284 -h -q 20 -L {input.bed} -o {output.bam} {input.bam}
        samtools index {output.bam}
        """
rule cfdna_wgs_samtools_alignment_qc:
    conda:
        "../config/cfdna-wgs-conda-env.yaml",
    input:
        f"{config['cfdna-wgs-dir']}/bams/{{library_id}}.{{ref_name}}.{{align_method}}.{{processing}}.bam",
    log:
        flagstat = f"{config['log-dir']}/{{library_id}}.{{ref_name}}.{{align_method}}.{{processing}}_cfdna_wgs_flagstat.log",
        samstat = f"{config['log-dir']}/{{library_id}}.{{ref_name}}.{{align_method}}.{{processing}}_cfdna_wgs_samstat.log",
    output:
        flagstat = f"{config['cfdna-wgs-dir']}/qc/{{library_id}}.{{ref_name}}.{{align_method}}.{{processing}}_flagstat.txt",
        samstat = f"{config['cfdna-wgs-dir']}/qc/{{library_id}}.{{ref_name}}.{{align_method}}.{{processing}}_samstat.txt",
    params:
        script = f"{config['cfdna-scripts-dir']}/samtools_alignment_qc.sh",
        threads = 8,
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
rule cfdna_wgs_mosdepth:
    conda:
        "../config/cfdna-wgs/conda-env.yaml",
    input:
        bam = f"{config['cfdna-wgs-dir']}/bams/{{library_id}}.{{ref_name}}.{{align_method}}.dedup.coorsort.filt.bam",
        index = f"{config['cfdna-wgs-dir']}/bams/{{library_id}}.{{ref_name}}.{{align_method}}.dedup.coorsort.filt.bam",
    output:
        summary = f"{config['cfdna-wgs-dir']}/qc/mosdepth_{{library_id}}.{{ref_name}}.{{align_method}}.mosdepth.summary.txt",
        global_dist = f"{config['cfdna-wgs-dir']}/qc/mosdepth_{{library_id}}.{{ref_name}}.{{align_method}}.mosdepth.global.dist.txt",
        region_dist = f"{config['cfdna-wgs-dir']}/qc/mosdepth_{{library_id}}.{{ref_name}}.{{align_method}}.mosdepth.region.dist.txt",
        regions = f"{config['cfdna-wgs-dir']}/qc/mosdepth_{{library_id}}.{{ref_name}}.{{align_method}}.regions.bed.gz",
        regions_idx = f"{config['cfdna-wgs-dir']}/qc/mosdepth_{{library_id}}.{{ref_name}}.{{align_method}}.regions.bed.gz.csi",
        quantized = f"{config['cfdna-wgs-dir']}/qc/mosdepth_{{library_id}}.{{ref_name}}.{{align_method}}.quantized.bed.gz",
        quantized_idx = f"{config['cfdna-wgs-dir']}/qc/mosdepth_{{library_id}}.{{ref_name}}.{{align_method}}.quantized.bed.gz.csi",
        thresholds = f"{config['cfdna-wgs-dir']}/qc/mosdepth_{{library_id}}.{{ref_name}}.{{align_method}}.thresholds.bed.gz",
        thresholds_idx = f"{config['cfdna-wgs-dir']}/qc/mosdepth_{{library_id}}.{{ref_name}}.{{align_method}}.thresholds.bed.gz.csi",
    params:
        script = f"{config['cfdna-scripts-dir']}/cdfna_wgs_mosdepth.sh",
        quant_levels = config["mosdepth-quant-levels"],
        out_dir = f"{config['cfdna-wgs-dir']}/qc",
    threads: 8,
    resources:
        concurrency = 20,
    shell:
        """
        {params.script} \
        {input.bam} \
        {params.out_dir} \
        {wildcards.library_id}.{wildcards.ref_name}.{wildcards.align_method} \
        '{params.quant_levels}' \
        {threads}
        """
# Get fragment sizes using deepTools
rule cfdna_wgs_frag_bampefragsize:
    conda:
        "../config/cfdna-wgs-conda-env.yaml",
    input:
        lambda wildcards: expand(f"{config['data-dir']}/cfdna_wgs/bams/{{library}}.{{ref_name}}.{{align_method}}.dedup.coorsort.filt.bam",
                                 library = cdfna_wgs_map[wildcards.lib_set]['libs'],
                                 build = lib_map[wildcards.lib_set]['build']),
    log: f"{config['log-dir']}/{{lib_set}}_cfdna_wgs_bampefragsize.log",
    output:
        raw = f"{config['data-dir']}/cfdna_wgs/qc/{{lib_set}}_bampefragsize.txt",
        hist = f"{config['data-dir']}/cfdna_wgs/qc/{{lib_set}}_bampefragsize.png",
    params:
        blacklist = lambda wildcards: lib_map[wildcards.lib_set]['blacklist'],
        script = f"{config['cfdna-scripts-dir']}/bampefragsize.sh",
    shell:
        """
        {params.script} \
        "{input}" \
        {log} \
        {output.hist} \
        {output.raw} \
        {params.blacklist} \
        12
        """
