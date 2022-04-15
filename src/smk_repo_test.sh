eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"

conda activate snakemake

output_dirs=( "bam" "processed-fastq" "qc" "symlink-fastq" "unpaired-fastq" )

for dir in ${output_dirs[@]};
do
               if [ -d test/${dir} ]; then \rm -rf test/${dir}; fi
done

snakemake \
    --configfile config/repo_test.yaml \
    --cores $threads \
    --rerun-incomplete \
    --use-singularity \
    --forceall \
    --snakefile ./workflow/read_preprocess.smk
