# Check for parameters, return usage if empty
if [[ $# -eq 0 ]] || [[ smk_dry_run.sh == "h" ]] ; then
    printf "\n usage: smk_forced_run.sh config_file smk_file
           \n "
else
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"

conda activate snakemake

snakemake \
    --configfile $1 --cores 4 \
    --use-singularity \
    --forceall \
    --printshellcmds \
    --snakefile $2
fi
