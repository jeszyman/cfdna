#!/bin/bash

eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"

echo "The following `*.smk` archives were found; select one:"

# set the prompt used by select, replacing "#?"
PS3="Use number to select an option"

select filename in ./workflow/*.smk

do
    if [[ "$filename" == "" ]]
    then
        echo "'$REPLY' is not a valid number"
        continue
    fi
    conda activate snakemake
    snakemake --dry-run --snakefile $filename \
              --configfile config/${HOSTNAME}.yaml \
              --cores $threads \
              --directory /drive3/users/jszymanski/repos/mpnst-preprocessing \
              --rerun-incomplete \
              --singularity-args "--bind $mntpt:$mntpt" \
              --use-singularity
    break
done
