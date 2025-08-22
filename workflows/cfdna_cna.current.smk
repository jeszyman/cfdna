rule fragment_filter:
    conda: f"{config['envs']['cfdna-cna']}"
    input:
        bam = f"{D_CFDNA_CNA}/bams/input/{{lib}}.input.bam"
    output:
        bam = f"{D_CFDNA_CNA}/bams/filt/{{lib}}.frag{{frag}}.bam",
        bai = f"{D_CFDNA_CNA}/bams/filt/{{lib}}.frag{{frag}}.bam.bai"
    log:
        f"{D_LOGS}/{{lib}}.frag{{frag}}_filter.log"
    params:
        script = f"{R_CFDNA}/scripts/cfdna_frag_filt.sh",
        tmp = f"{D_LOCAL_TMP}/frag"
    threads: 12
    resources:
        concurrency = 25
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.bam})" "{params.tmp}/{wildcards.lib}"
        if [ "{wildcards.frag}" = "_all" ]; then
            cp {input.bam} {output.bam}
            samtools index -@ {threads} {output.bam}
        else
            {params.script} \
              {input.bam} \
              {params.tmp}/{wildcards.lib}.frag{wildcards.frag}.nohead \
              $(echo {wildcards.frag} | cut -d_ -f1) \
              $(echo {wildcards.frag} | cut -d_ -f2) \
              {threads} \
              {params.tmp}/{wildcards.lib}/frag{wildcards.frag}.onlyhead \
              {output.bam} \
              {params.tmp} &> {log}
            samtools index -@ {threads} {output.bam}
        fi
        """
rule cfdna_cna_downsample_bam:
    conda: f"{config['envs']['cfdna-cna']}"
    input:
        bam = f"{D_CFDNA_CNA}/bams/filt/{{lib}}.frag{{frag}}.bam",
        bai = f"{D_CFDNA_CNA}/bams/filt/{{lib}}.frag{{frag}}.bam.bai",
    log:
        f"{D_LOGS}/{{lib}}.frag{{frag}}.ds{{mil_reads}}_cfdna_cna_downsample_bam.log",
    output:
        bam = f"{D_CFDNA_CNA}/bams/ds/{{lib}}.frag{{frag}}.ds{{mil_reads}}.bam",
        bai = f"{D_CFDNA_CNA}/bams/ds/{{lib}}.frag{{frag}}.ds{{mil_reads}}.bam.bai",
    params:
        milreads = lambda wildcards: wildcards.mil_reads,
        script = f"{config['cfdna-scripts-dir']}/downsample_bam.sh",
    shell:
        r"""
        if [ "{wildcards.mil_reads}" = "_none" ]; then
            cp {input.bam} {output.bam}
            samtools index {output.bam}
        else
            {params.script} \
              {input.bam} \
              {params.milreads} \
              {output.bam} &> {log}
            samtools index {output.bam}
        fi
        """
# WIG generation - handles all BAM types
rule make_wig:
    conda: f"{config['envs']['cfdna-cna']}"
    input:
        f"{D_CFDNA_CNA}/bams/ds/{{lib}}.frag{{frag}}.ds{{mil_reads}}.bam",
    output:
        f"{D_CFDNA_CNA}/wigs/{{lib}}.frag{{frag}}.ds{{mil_reads}}.wig",
    log:
        f"{D_LOGS}/{{lib}}.frag{{frag}}.ds{{mil_reads}}_make_wig.log",
    params:
        window = "1000000",
        quality = 20,
        ichor_wig_dir = f"{D_CFDNA_CNA}/wigs",
    shell:
        """
        mkdir -p "{params.ichor_wig_dir}"
        readCounter \
        --window {params.window} \
        --quality {params.quality} \
	--chromosome "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY" \
        {input} > {output}
        """
rule cfdna_ichor_pon_list:
    conda: f"{config['envs']['cfdna-cna']}"
    input:
        wigs=lambda wc: [
            f"{D_CFDNA_CNA}/wigs/{lib}.frag{PON_SETS[wc.pon_id]['frag']}.ds{PON_SETS[wc.pon_id]['ds']}.wig"
            for lib in PON_SETS[wc.pon_id]["libs"]
        ]
    output:
        f"{D_CFDNA_CNA}/pon/{{pon_id}}_normals_list.txt"
    log:
        f"{D_LOGS}/pon_{{pon_id}}_wiglist.log"
    shell:
        r"""
        rm -f {output}
        for f in {input}; do
          [ -s "$f" ] && readlink -f "$f"
        done | sort -u > {output}
        if [ ! -s {output} ]; then
          echo "No WIGs found for pon_id={wildcards.pon_id}" >&2
          exit 2
        fi
        """

rule cfdna_ichor_pon:
    conda: f"{config['envs']['cfdna-cna']}"
    input:
        f"{D_CFDNA_CNA}/pon/{{pon_id}}_normals_list.txt",
    output:
        txt = f"{D_CFDNA_CNA}/pon/{{pon_id}}_median.txt",
        rds = f"{D_CFDNA_CNA}/pon/{{pon_id}}_median.rds",
    log:
        f"{D_LOGS}/pon_{{pon_id}}_build.log",
    params:
        out_prefix = lambda wc: f"{D_CFDNA_CNA}/pon/{wc.pon_id}",
        ichor_repo = config["ichor_repo"],
    shell:
        r"""
        set -euo pipefail
        Rscript {params.ichor_repo}/scripts/createPanelOfNormals.R \
          --filelist {input} \
          --chrNormalize "paste0('chr',1:22)" \
          --chrs "c(paste0('chr',1:22),'chrX')" \
          --gcWig {params.ichor_repo}/inst/extdata/gc_hg38_1000kb.wig \
          --mapWig {params.ichor_repo}/inst/extdata/map_hg38_1000kb.wig \
          --centromere {params.ichor_repo}/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt \
          --libdir {params.ichor_repo} \
          --outfile {params.out_prefix} &> {log}
        # Ensure both artifacts exist (script should produce both)
        test -s {output.txt} && test -s {output.rds}
        """
rule cfdna_ichor:
    conda: f"{config['envs']['cfdna-cna']}"
    input:
        wig = f"{D_CFDNA_CNA}/wigs/{{lib}}.frag{{frag}}.ds{{mil_reads}}.wig",
    output:
        seg = f"{D_CFDNA_CNA}/ichor/{{ichor_set}}/{{lib}}.frag{{frag}}.ds{{mil_reads}}.ichor_{{ichor_set}}.cna.seg",
    params:
        ichor_out_main_dir = lambda wc: f"{D_CFDNA_CNA}/ichor/{wc.ichor_set}",
        ichor_repo   = config["ichor_repo"],
        genome       = lambda wc: ichor_map[wc.ichor_set]['genome'],
        includeHOMD  = lambda wc: ichor_map[wc.ichor_set]['includeHOMD'],
        maxCN        = lambda wc: ichor_map[wc.ichor_set]['maxCN'],
        normal_prior = lambda wc: ichor_map[wc.ichor_set]['normal_prior'],
        ploidy_prior = lambda wc: ichor_map[wc.ichor_set]['ploidy_prior'],
        chrs         = lambda wc: ichor_map[wc.ichor_set]['chrs'],
        chrTrain     = lambda wc: ichor_map[wc.ichor_set]['chrTrain'],
        chrNormalize = lambda wc: ichor_map[wc.ichor_set]['chrNormalize'],
        txnE         = lambda wc: ichor_map[wc.ichor_set]['txnE'],
        txnStrength  = lambda wc: ichor_map[wc.ichor_set]['txnStrength'],
        estimateNormal        = lambda wc: ichor_map[wc.ichor_set]['estimateNormal'],
        estimatePloidy        = lambda wc: ichor_map[wc.ichor_set]['estimatePloidy'],
        estimateScPrevalence  = lambda wc: ichor_map[wc.ichor_set]['estimateScPrevalence'],
        normal_panel_opt = lambda wc: (
            (lambda repo, val: (
                "" if (val is None or val == "" or val == "OFF") else
                f'--normalPanel "{repo}/inst/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds"' if val == "DEFAULT" else
                (f'--normalPanel "{val}"' if val.startswith(("/", "~")) else
                 f'--normalPanel "{repo}/{val}"')
            ))(config.get("ichor_repo", config.get("ichor-repo")), ichor_map[wc.ichor_set].get("normalPanel"))
        ),
    log:
        f"{D_LOGS}/{{lib}}.frag{{frag}}.ds{{mil_reads}}.{{ichor_set}}_ichor.log",
    shell:
        r"""
        mkdir -p {params.ichor_out_main_dir}
        Rscript {params.ichor_repo}/scripts/runIchorCNA.R \
          --id {wildcards.lib}.frag{wildcards.frag}.ds{wildcards.mil_reads}.ichor_{wildcards.ichor_set} \
          --WIG {input.wig} \
          --genomeBuild "{params.genome}" \
          --gcWig {params.ichor_repo}/inst/extdata/gc_hg38_1000kb.wig \
          --mapWig {params.ichor_repo}/inst/extdata/map_hg38_1000kb.wig \
          --centromere {params.ichor_repo}/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt \
          {params.normal_panel_opt} \
          --includeHOMD {params.includeHOMD} \
          --ploidy "{params.ploidy_prior}" \
          --normal "{params.normal_prior}" \
          --maxCN {params.maxCN} \
          --chrs "{params.chrs}" \
          --chrTrain "{params.chrTrain}" \
          --chrNormalize "{params.chrNormalize}" \
          --estimateNormal {params.estimateNormal} \
          --estimatePloidy {params.estimatePloidy} \
          --estimateScPrevalence {params.estimateScPrevalence} \
          --txnE {params.txnE} \
          --txnStrength {params.txnStrength} \
          --outDir {params.ichor_out_main_dir} \
          --libdir {params.ichor_repo}
        """
