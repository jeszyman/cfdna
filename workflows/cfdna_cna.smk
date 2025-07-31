rule cfdna_frag_filt:
    conda:
        f"{config['cfdna-cna-conda-env']}",
    input:
        f"{config['cfdna-cna-dir']}/input_bams/{{wkflow_id}}.bam",
    log:
        f"{config['log-dir']}/{{wkflow_id}}_{{frag_distro}}_cfdna_frag_filt.log",
    output:
        nohead = temp(f"{config['cfdna-cna-dir']}/bams/{{wkflow_id}}.frag{{frag_distro}}.nohead"),
        onlyhead = temp(f"{config['cfdna-cna-dir']}/bams/{{wkflow_id}}.frag{{frag_distro}}.onlyhead"),
        final = f"{config['cfdna-cna-dir']}/bams/{{wkflow_id}}.frag{{frag_distro}}.bam",
        index = f"{config['cfdna-cna-dir']}/bams/{{wkflow_id}}.frag{{frag_distro}}.bam.bai",
    params:
        script = f"{config['cfdna-scripts-dir']}/cfdna_frag_filt.sh",
        tmp_dir = config["local-tmp-dir"],
    threads:
        12
    shell:
        """
        /home/ext_szymanski_jeffrey_mayo_edu/repos/cfdna/scripts/cfdna_frag_filt.sh \
        {input} \
        {output.nohead} \
        $(echo {wildcards.frag_distro} | cut -d_ -f1) \
        $(echo {wildcards.frag_distro} | cut -d_ -f2) \
        {threads} \
        {output.onlyhead} \
        {output.final} \
        {params.tmp_dir} &> {log}
        samtools index {output.final}
        """
rule cfdna_cna_downsample_bam:
    conda:
        f"{config['cfdna-cna-conda-env']}",
    input:
        f"{config['cfdna-cna-dir']}/bams/{{wkflow_id}}.frag{{frag_distro}}.bam",
    log:
        f"{config['log-dir']}/{{wkflow_id}}.frag{{frag_distro}}_{{mil_reads}}_cfdna_cna_downsample.log",
    output:
        bam = f"{config['cfdna-cna-dir']}/ds-bams/{{wkflow_id}}.frag{{frag_distro}}.ds{{mil_reads}}.bam",
        bai = f"{config['cfdna-cna-dir']}/ds-bams/{{wkflow_id}}.frag{{frag_distro}}.ds{{mil_reads}}.bam.bai",
    params:
        milreads = lambda wildcards: wildcards.mil_reads,
        script = f"{config['cfdna-scripts-dir']}/downsample_bam.sh",
    shell:
        """
        {params.script} \
        {input} \
        {params.milreads} \
        {output.bam} &> {log}
        """
rule make_wig:
    conda:
        f"{config['cfdna-cna-conda-env']}",
    input:
        bam = f"{config['cfdna-cna-dir']}/ds-bams/{{wkflow_id}}.frag{{frag_distro}}.ds{{mil_reads}}.bam",
        bai = f"{config['cfdna-cna-dir']}/ds-bams/{{wkflow_id}}.frag{{frag_distro}}.ds{{mil_reads}}.bam.bai",
    log:
        f"{config['log-dir']}/{{wkflow_id}}.frag{{frag_distro}}.ds{{mil_reads}}_cfdna_cna_make_wig.log",
    output:
        wig = f"{config['cfdna-cna-dir']}/wigs/{{wkflow_id}}.frag{{frag_distro}}.ds{{mil_reads}}.wig",
    params:
        window = "1000000",
        quality = 20,
        ichor_wig_dir = f"{config['cfdna-cna-dir']}/wigs",
    shell:
        """
        mkdir -p "{params.ichor_wig_dir}"
        readCounter \
        --window {params.window} \
        --quality {params.quality} \
	--chromosome "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY" \
        {input} > {output}
        """
rule cfdna_ichor:
    conda:
        f"{config['cfdna-cna-conda-env']}",
    input:
        f"{config['cfdna-cna-dir']}/wigs/{{wkflow_id}}.frag{{frag_distro}}.ds{{mil_reads}}.wig",
    output:
        f"{config['cfdna-cna-dir']}/ichor/{{wkflow_id}}.frag{{frag_distro}}.ds{{mil_reads}}/{{wkflow_id}}.frag{{frag_distro}}.ds{{mil_reads}}.cna.seg",
    params:
        ichor_out_main_dir = f"{config['cfdna-cna-dir']}/ichor",
        ichor_repo = config["ichor_repo"],
    shell:
        """
        mkdir -p $(dirname {params.ichor_out_main_dir}/{wildcards.wkflow_id}.frag{wildcards.frag_distro}.ds{wildcards.mil_reads})
        Rscript {params.ichor_repo}/scripts/runIchorCNA.R \
            --id {wildcards.wkflow_id}.frag{wildcards.frag_distro}.ds{wildcards.mil_reads} \
            --WIG {input} \
            --normal "c(0.95, 0.99, 0.995, 0.999)" \
            --genomeBuild hg38 \
            --ploidy "c(2)" \
            --gcWig {params.ichor_repo}/inst/extdata/gc_hg38_1000kb.wig \
            --mapWig {params.ichor_repo}/inst/extdata/map_hg38_1000kb.wig \
            --centromere {params.ichor_repo}/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt \
            --normalPanel {params.ichor_repo}/inst/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds \
            --includeHOMD FALSE \
            --chrs "c(1:22)" \
            --chrTrain "c(1:22)" \
            --estimateNormal TRUE \
            --estimatePloidy TRUE \
            --estimateScPrevalence TRUE \
            --scStates "c()" \
            --txnE 0.9999 \
            --txnStrength 10000 \
            --outDir {params.ichor_out_main_dir}/{wildcards.wkflow_id}.frag{wildcards.frag_distro}.ds{wildcards.mil_reads} \
            --libdir {params.ichor_repo}
        """
# Snakefile
# at the top of your Snakefile (if not already):
from pathlib import Path
import os

# … your other rules …

rule cfdna_cna_extract_tumor_fractions:
    conda:
        f"{config['cfdna-cna-conda-env']}"
    input:
        # dynamically find every .params.txt under ichor/
        params=lambda wc: sorted(
            str(p) for p in Path(config['cfdna-cna-dir'], "ichor")
                       .rglob("*.params.txt")
        )
    output:
        tf_summary=os.path.join(
            config['cfdna-cna-dir'], "ichor", "ichor_tumor_fractions.tsv"
        )
    shell:
        r"""
        # write header
        echo -e "library\ttf" > {output.tf_summary}

        # for each params.txt, pull the second field on line 2
        for f in {input.params}; do
            sample=$(basename "$f" .params.txt)
            tf=$(awk 'NR==2 {{print $2}}' "$f")
            echo -e "$sample\t$tf"
        done >> {output.tf_summary}
        """
