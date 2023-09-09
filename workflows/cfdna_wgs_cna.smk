#########1#########2#########3#########4#########5#########6#########7#########8
#                                                                              #
#                   Snakefile for Analysis of Cell-free DNA                    #
#    Whole Genome Sequencing Copy Number Alteration and Fragmentomics          #
#                                                                              #
#########1#########2#########3#########4#########5#########6#########7#########8

rule cfdna_wgs_downsample_bam:
    input: f"{cfdna_wgs_dir}/bams/{{library}}_{{build}}_filt.bam",
    log: f"{log_dir}/{{library}}_{{build}}_ds{{mil_reads}}_cfdna_wgs_downsample_bam.log",
    output: f"{cfdna_wgs_dir}/bams/{{library}}_{{build}}_ds{{mil_reads}}.bam",
    params:
        milreads = lambda wildcards: wildcards.mil_reads,
        script = f"{cfdna_script_dir}/downsample_bam.sh",
        threads = threads,
    shell:
        """
        {params.script} \
        {input} \
        {params.milreads} &> {log}
        """

rule cfdna_wgs_frag_filt:
    input: f"{cfdna_wgs_dir}/bams/{{library}}_{{build}}_ds{{mil_reads}}.bam",
    log: f"{log_dir}/{{library}}_{{build}}_ds{{mil_reads}}_{{frag_distro}}_cfdna_wgs_frag_filt.log",
    params: script = f"{cfdna_script_dir}/wgs_frag_filt.sh",
    output:
        nohead = temp(f"{cfdna_wgs_dir}/frag_bams/{{library}}_{{build}}_ds{{mil_reads}}_frag{{frag_distro}}.nohead"),
        onlyhead = temp(f"{cfdna_wgs_dir}/frag_bams/{{library}}_{{build}}_ds{{mil_reads}}_frag{{frag_distro}}.onlyhead"),
        final = f"{cfdna_wgs_dir}/frag_bams/{{library}}_{{build}}_ds{{mil_reads}}_frag{{frag_distro}}.bam",
        index = f"{cfdna_wgs_dir}/frag_bams/{{library}}_{{build}}_ds{{mil_reads}}_frag{{frag_distro}}.bam.bai",
    shell:
        """
        frag_min=$(echo {wildcards.frag_distro} | sed -e "s/_.*$//g")
        frag_max=$(echo {wildcards.frag_distro} | sed -e "s/^.*_//g")
        {params.script} \
        {input} \
        {output.nohead} \
        $frag_min \
        $frag_max \
        4 \
        {output.onlyhead} \
        {output.final} &> {log}
        samtools index {output.final}
        """

rule bam_to_wig:
    input: f"{cfdna_wgs_dir}/frag_bams/{{library}}_{{build}}_ds{{mil_reads}}_frag{{frag_distro}}.bam",
    log: f"{log_dir}/{{library}}_{{build}}_{{mil_reads}}_{{frag_distro}}_cfdna_bam_to_wig.log",
    output: ensure(f"{cfdna_wgs_dir}/cna/wigs/{{library}}_{{build}}_ds{{mil_reads}}_frag{{frag_distro}}.wig", non_empty=True),
    params:
        chrs = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY",
        out_dir = f"{cfdna_wgs_dir}/cna/wigs",
    shell:
        """
        mkdir -p {params.out_dir} && readCounter --window 1000000 --quality 20 --chromosome {params.chrs} {input} > {output}
        """

rule cfdna_wgs_ichor:
    input: f"{cfdna_wgs_dir}/cna/wigs/{{library}}_{{build}}_ds{{mil_reads}}_frag{{frag_distro}}.wig",
    log: f"{log_dir}/{{library}}_{{build}}_ds{{mil_reads}}_frag{{frag_distro}}_cfdna_wgs_ichor.log",
    output: f"{cfdna_wgs_dir}/cna/ichor_nopon/{{library}}_{{build}}_ds{{mil_reads}}_frag{{frag_distro}}.RData",
    params: out_dir = f"{cfdna_wgs_dir}/cna/ichor_nopon",
    shell:
        """
        Rscript /opt/miniconda3/envs/mpnst/bin/runIchorCNA.R \
        --id {wildcards.library}_frag{wildcards.frag_distro} \
        --WIG {input} \
        --gcWig /opt/miniconda3/envs/mpnst/lib/R/library/ichorCNA/extdata/gc_hg38_1000kb.wig \
        --mapWig /opt/miniconda3/envs/mpnst/lib/R/library/ichorCNA/extdata/map_hg38_1000kb.wig \
        --centromere /opt/miniconda3/envs/mpnst/lib/R/library/ichorCNA/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt \
        --normal "c(0.95, 0.99, 0.995, 0.999)" \
        --ploidy "c(2)" \
        --maxCN 3 \
        --estimateScPrevalence FALSE \
        --scStates "c()" \
        --outDir {params.out_dir}
        """
