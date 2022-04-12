eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"

conda activate snakemake
source config/${HOSTNAME}.sh

snakemake \
    --configfile config/repo_test.yaml \
    --cores $threads \
    --directory ${repo} \
    --dry-run \
    --rerun-incomplete \
    --use-singularity \
    --snakefile ./workflow/read_qc.smk

    &&
    
    snakemake \
        --configfile config/repo_test.yaml \
        --cores $threads \
        --directory ${repo} \
        --rerun-incomplete \
        --use-singularity \
        -F \
        --snakefile ./workflow/read_qc.smk
