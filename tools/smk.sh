#!/usr/bin/env bash

#######################################
###   Snakemake Conda Run Wrapper   ###
#######################################

#bash.usage

# Variables
env="${1}"
config="${2}"
snkfile="${3}"

# Necessary to run conda snakemake command in shell script
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"

conda activate $env

run_dry(){
    snakemake \
        --configfile $config \
        --cores 1 \
        --dry-run \
        --printshellcmds \
        --snakefile $snkfile
}

run(){
    snakemake \
        --configfile $config \
        --cores $nproc \
        --printshellcmds \
        --snakefile $snkfile
}

if [ "$4" == "run" ]; then
    run $env $config $snkfile
else
    run_dry $env $config $snkfile
fi
