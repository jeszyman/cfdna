# [[file:~/repos/cfdna-wgs/cfdna-wgs.org::*Repository%20setup%20and%20administration][Repository setup and administration:3]]
#!/usr/bin/env bash

print_usage(){
    cat <<- EOF

usage: conda_yaml_build_and_update.sh <CONDA ENV YAML FILE>

Conda environment builder and updater
Assumes mamba is present in your base conda environment.

example: conda_yaml_build_and_update.sh test_env.yaml

EOF
}

# Path to the YAML file
yaml_path=$1

# Necessary to run conda in shell script
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"

conda activate base

build_or_update(){
    # Name of the Conda environment
    env_name=$(head $yaml_path -n 1 | sed 's/^.*:.//g')
    # Check if the environment already exists
    if conda env list | grep -q "^$env_name\s"; then
        # Update the existing environment
        mamba env update -n "$env_name" -f "$yaml_path"
    else
        # Create a new environment
        mamba env create -n "$env_name" -f "$yaml_path"
    fi
}

conda deactivate

if [ $# -ne 1 ]; then print_usage && exit 1; fi

build_or_update
# Repository setup and administration:3 ends here
