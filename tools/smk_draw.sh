eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"

conda activate snakemake

snakemake \
    --configfile config/repo_test.yaml \
    --cores $threads \
    --rulegraph \
    --snakefile ./workflow/read_preprocess.smk | dot -Tpdf > resources/read_preprocess_dagtmp/test.pdf
