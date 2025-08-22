import pandas as pd
import os

configfile: "./config/TEST_cfdna_cna.yaml"

shell.executable("/bin/bash")

def resolve_config_paths(config_dict):
    for k, v in config_dict.items():
        if isinstance(v, str):
            config_dict[k] = os.path.expandvars(os.path.expanduser(v))
        elif isinstance(v, dict):
            resolve_config_paths(v)
        elif isinstance(v, list):
            config_dict[k] = [os.path.expandvars(os.path.expanduser(i)) if isinstance(i, str) else i for i in v]

resolve_config_paths(config)

# Directories
D_PROJECT = config['project-data-dir']
D_CFDNA_CNA = f"{D_PROJECT}/cfdna-cna"
D_LOGS = f"{D_PROJECT}/logs"
D_LOCAL_TMP = f"{D_PROJECT}/tmp"
R_CFDNA = config['cfdna-repo']

# Sample info
SAMPLES = config["samples-tsv"]
LIBS = sorted(pd.read_csv(SAMPLES, sep="\t")["library_id"].unique())

# Parameters
FRAGS = config.get("frag_distros", ["90_150", "_all"])
MILS = config.get("ds_mil_reads", ["2", "_none"])

# ichor_map presets consumed by the Snakemake rule via {wildcards.ichor_set}
# Values are strings when interpolated directly into the CLI.
ichor_map = {
    'base': {
        # Genome / scope
        'genome': "hg38",                       # genomeBuild passed to ichorCNA

        # Copy-state space
        'includeHOMD': "FALSE",                 # omit CN=0 state (reduces spurious deep losses at low TF)
        'maxCN': 4,                             # cfDNA-friendly cap; limits overfitting / WGD “explanations”

        # Priors (fit TF around high-normal; keep baseline near diploid)
        'normal_prior': "c(0.60,0.70,0.80,0.90,0.95,0.97,0.99)",  # cheap to try wide; pick best by likelihood
        'ploidy_prior': "c(2)",                 # search around 2N; prevents drift toward ~3N

        # Chromosome sets
        'chrs': "c(paste0('chr',1:22),'chrX')", # segment/report X, but…
        'chrTrain': "paste0('chr',1:22)",       # …train HMM on autosomes only (stable baseline)
        'chrNormalize': "paste0('chr',1:22)",   # …normalize on autosomes only (avoid sex-chrom bias)

        # Panel of Normals
        'normalPanel': None,                    # None → omit --normalPanel (PoN OFF)

        # HMM smoothing
        'txnE': "0.995",                        # transition prob (higher = smoother, fewer breakpoints)
        'txnStrength': "250",                   # segmentation penalty (higher = smoother)

        # Estimation flags
        'estimateNormal': "TRUE",               # fit normal fraction (i.e., TF)
        'estimatePloidy': "TRUE",               # allow small drift around 2 (e.g., 2.02)
        'estimateScPrevalence': "FALSE",        # keep subclones off until clonal fit is clean
    },

    'cfDNA_PoN_guard': {
        # Same as base but with a *guarded* PoN + tighter normal grid
        'genome': "hg38",

        'includeHOMD': "FALSE",
        'maxCN': 4,

        'normal_prior': "c(0.90,0.95,0.97,0.99)",  # high-normal grid stabilizes TF at low fractions
        'ploidy_prior': "c(2)",

        'chrs': "c(paste0('chr',1:22),'chrX')",
        'chrTrain': "paste0('chr',1:22)",
        'chrNormalize': "paste0('chr',1:22)",

        'normalPanel': 'inst/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds',
        # NOTE: Use this only if PoN is processed identically to samples (same bins, fragment window, chr style).

        'txnE': "0.995",
        'txnStrength': "250",

        'estimateNormal': "TRUE",
        'estimatePloidy': "TRUE",
        'estimateScPrevalence': "FALSE",
    },

    'ha_low_tumor': {

        # Genome / naming (stock ichor defaults)
        'genome': "hg38",

        # Copy-state space
        'includeHOMD': "FALSE",           # default
        'maxCN': 3,

        # Priors
        'normal_prior': "c(0.5)",         # default (single start at 50% normal)
        'ploidy_prior': "c(0.95, 0.99, 0.995, 0.999)",

        # Chromosome sets (stock uses non-‘chr’ names)
        'chrs': "c(paste0('chr',1:22),'chrX')",
        'chrTrain': "paste0('chr',1:22)",
        'chrNormalize': "paste0('chr',1:22)",

        # Panel of Normals
        'normalPanel': "DEFAULT",

        # HMM smoothing (very aggressive by default)
        'txnE': "0.9999999",              # default
        'txnStrength': "10000000",        # default (1e7)

        # Estimation flags
        'estimateNormal': "TRUE",         # default
        'estimatePloidy': "TRUE",         # default
        'estimateScPrevalence': "FALSE",
    },
        'TEST': {

        # Genome / naming (stock ichor defaults)
        'genome': "hg38",

        # Copy-state space
        'includeHOMD': "FALSE",           # default
        'maxCN': 3,

        # Priors
        'normal_prior': "c(0.5)",         # default (single start at 50% normal)
        'ploidy_prior': "c(0.95, 0.99, 0.995, 0.999)",

        # Chromosome sets (stock uses non-‘chr’ names)
        'chrs': "c(paste0('chr',1:22),'chrX')",
        'chrTrain': "paste0('chr',1:22)",
        'chrNormalize': "paste0('chr',1:22)",

        # Panel of Normals
        'normalPanel': f"{D_CFDNA_CNA}/pon/test_median.rds",

        # HMM smoothing (very aggressive by default)
        'txnE': "0.9999999",              # default
        'txnStrength': "10000000",        # default (1e7)

        # Estimation flags
        'estimateNormal': "TRUE",         # default
        'estimatePloidy': "TRUE",         # default
        'estimateScPrevalence': "FALSE",
    },

}



# --- Preamble: define PoN sets here (no config) ------------------------------
# Each key is a pon_id that will produce D_CFDNA_CNA/pon/{pon_id}_median.txt
PON_SETS = {
    "test": {
        "libs": ["kh_01", "kh_02"],  # library_id values
        "frag": "90_150",
        "ds": "_none",
    },
    # add more sets if needed
    # "pon2": {"libs": [...], "frag": "90_150", "ds": "_none"},
}



wildcard_constraints:
    frag=r"(?:_all|\d+_\d+)",
    mil_reads=r"(?:_none|\d+)"

# All final outputs
rule all:
    input:
        # # Fragment filter
        # expand(f"{D_CFDNA_CNA}/bams/filt/{{lib}}.frag{{frag}}.bam",
        #        lib=LIBS, frag=FRAGS),
        # #
        # # Downsample
        # expand(f"{D_CFDNA_CNA}/bams/ds/{{lib}}.frag{{frag}}.ds{{mil_reads}}.bam",
        #        lib=LIBS, frag=FRAGS, mil_reads=MILS),
        #
        # Bam -> Wiggle
        # expand(f"{D_CFDNA_CNA}/wigs/{{lib}}.frag{{frag}}.ds{{mil_reads}}.wig",
        #        lib=LIBS, frag=FRAGS, mil_reads=MILS),
        #
        # Make PoNs
        expand(f"{D_CFDNA_CNA}/pon/{{pon_id}}_normals_list.txt",
               pon_id=["test"]),
        expand(f"{D_CFDNA_CNA}/pon/{{pon_id}}_median.rds",
               pon_id=["test"]),

        expand(f"{D_CFDNA_CNA}/pon/{{pon_id}}_median.txt",
               pon_id=["test"]),

        #
        # NEW ichor
        expand(f"{D_CFDNA_CNA}/ichor/{{ichor_set}}/{{lib}}.frag{{frag}}.ds{{mil_reads}}.ichor_{{ichor_set}}.cna.seg",
               lib=LIBS, frag=FRAGS, mil_reads=MILS, ichor_set=['base','cfDNA_PoN_guard','default','ha_low_tumor', 'TEST']),

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
