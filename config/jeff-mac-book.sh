#!/usr/bin/env bash
#########1#########2#########3#########4#########5#########6#########7#########8

#####################################################################
###   Bash configuration script for MPNST fragmentomics project   ###
#####################################################################

# Host-local variables
data_dir=/mnt/ris/aadel/mpnst
mntpt=/mnt/ris/aadel
repo=/home/jeszyman/repos/mpnst-preprocessing
sif_dir=/home/jeszyman/sing_containers
threads=8

# Check git file hook is read-able
if [ -r "${repo}/.git/hooks/precommit" ]; then
   echo "Git size check is read-able"
else
    echo
    "Git size check is not read-able"
    exit 1
fi
          
# Check mount point  
if grep -qs $mntpt /proc/mounts; then
    echo "RIS storage mounted."
else
    echo "RIS storage NOT mounted!!!"
fi

# Check local singularity build present
if [ -r ${sif_dir}/frag.sif ]; then
    echo "Local singularity file present"
else
    echo "Local singularity file not found."
fi

# Check if local singularity container is present and up-to-date
if [ -r "${HOME}/sing_containers/biotools.sif" ]; then
    echo "Local biotools singularity container exists"
else
    echo "No local biotools singularity container, attempting to fetch..."
    mkdir -p "${HOME}/sing_containers"
    cp /mnt/ris/aadel/jeszyman/sing_containers/biotools.sif "${HOME}/sing_containers"
fi 

if [ ${HOME}/sing_containers/biotools.sif -ot /mnt/ris/aadel/jeszyman/sing_containers/biotools.sif ];
then
    echo "Local biotools container out of date, updating..."
    cp /mnt/ris/aadel/jeszyman/sing_containers/biotools.sif "${HOME}/sing_containers"
else
    echo "Local biotools singularity container is up to date"
fi

sing_biotools() {
    singularit shell --bind /mnt:/mnt ~/sing_containers/biotools.sif            
}

launch_frag() { 
    if [ -f /.dockerenv ]; then
        echo "shell already in docker, exiting";
        exit 1;
    else
        docker run --env HOME=${HOME} --hostname ${HOSTNAME} --interactive --tty --volume /home/:/home/ --volume /tmp/:/tmp/ --volume /mnt/:/mnt/ --user $(id -u ${USER}) -w "$repo" jeszyman/frag /bin/bash;
    fi
}
