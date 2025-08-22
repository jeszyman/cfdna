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
        # Fragment filter
        expand(f"{D_CFDNA_CNA}/bams/filt/{{lib}}.frag{{frag}}.bam",
               lib=LIBS, frag=FRAGS),
        #
        # Downsample
        expand(f"{D_CFDNA_CNA}/bams/ds/{{lib}}.frag{{frag}}.ds{{mil_reads}}.bam",
               lib=LIBS, frag=FRAGS, mil_reads=MILS),
        #
        #Bam -> Wiggle
        expand(f"{D_CFDNA_CNA}/wigs/{{lib}}.frag{{frag}}.ds{{mil_reads}}.wig",
               lib=LIBS, frag=FRAGS, mil_reads=MILS),
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
               lib=LIBS, frag=FRAGS, mil_reads=MILS, ichor_set=['base','cfDNA_PoN_guard','ha_low_tumor', 'TEST']),


include: f"{config['cfdna-repo']}/workflows/cfdna_cna.current.smk"
