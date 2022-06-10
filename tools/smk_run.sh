#!/bin/bash
#########1#########2#########3#########4#########5#########6#########7#########8
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
####################################
###   Choose and Run Snakefile   ###
####################################

# Setup
#set -euxov pipefail

source config/${HOSTNAME}.sh
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
    echo $filename
    select run_option in dry_run normal force_final force_all
    do
        echo selected $run_option
        case $run_option in
            dry_run)
                conda activate snakemake
                snakemake \
                    --configfile config/${HOSTNAME}.yaml \
                    --cores $threads \
                    --directory ${repo} \
                    --dry-run \
                    --rerun-incomplete \
                    --singularity-args "--bind $mntpt:$mntpt" \
                    --use-singularity \
                    --snakefile $filename
                ;;
            normal) 
                conda activate snakemake
                select nohup_option in no yes
                do
                    case $nohup_option in
                        no)
                            snakemake \
                                --configfile config/${HOSTNAME}.yaml \
                                --cores $threads \
                                --directory ${repo} \
                                --singularity-args "--bind $mntpt:$mntpt" \
                                --use-singularity \
                                --snakefile $filename
                            ;;
                        yes)
                            nohup snakemake \
                                  --configfile config/${HOSTNAME}.yaml \
                                  --cores $threads \
                                  --directory ${repo} \
                                  --singularity-args "--bind $mntpt:$mntpt" \
                                  --use-singularity \
                                  --snakefile $filename
                            ;;
                    esac
                done                
                ;;
            force_final)
                conda activate snakemake
                select nohup_option in no yes
                do
                    case $nohup_option in
                        no)                                          
                            snakemake \
                                --configfile config/${HOSTNAME}.yaml \
                                --cores $threads \
                                --directory ${repo} \
                                --force \
                                --singularity-args "--bind $mntpt:$mntpt" \
                                --use-singularity \
                                --snakefile $filename
                            ;;
                        yes)
                            nohup snakemake \
                                  --configfile config/${HOSTNAME}.yaml \
                                  --cores $threads \
                                  --directory ${repo} \
                                  --force \
                                  --singularity-args "--bind $mntpt:$mntpt" \
                                  --use-singularity \
                                  --snakefile $filename
                            ;;
                    esac
                done                
                ;;            
            force_all)
                conda activate snakemake
                select nohup_option in no yes
                do
                    case $nohup_option in
                        no)                                          
                
                snakemake \
                    --configfile config/${HOSTNAME}.yaml \
                    --cores $threads \
                    --directory ${repo} \
                    -F \
                    --singularity-args "--bind $mntpt:$mntpt" \
                    --use-singularity \
                    --snakefile $filename
                                            ;;
                        yes)
                nohup snakemake \
                    --configfile config/${HOSTNAME}.yaml \
                    --cores $threads \
                    --directory ${repo} \
                    -F \
                    --singularity-args "--bind $mntpt:$mntpt" \
                    --use-singularity \
                    --snakefile $filename
                            ;;
                    esac
                done                
                ;;            
        esac
        break
    done
    break
done

#!/bin/bash
#########1#########2#########3#########4#########5#########6#########7#########8

####################################
###   Choose and Run Snakefile   ###
####################################

# Setup
#set -euxov pipefail
source config/${HOSTNAME}.sh
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
    echo $filename
    select run_option in dry_run normal force_final force_all
    do
        echo selected $run_option
        case $run_option in
            dry_run)
                source activate snakemake
                snakemake \
                    --configfile config/${HOSTNAME}.yaml \
                    --cores $threads \
                    --directory ${repo} \
                    --dry-run \
                    --rerun-incomplete \
                    --singularity-args "--bind $mntpt:$mntpt" \
                    --use-singularity \
                    --snakefile $filename
                ;;
            normal) 
                source activate snakemake
                select nohup_option in no yes
                do
                    case $nohup_option in
                        no)
                            snakemake \
                                --configfile config/${HOSTNAME}.yaml \
                                --cores $threads \
                                --directory ${repo} \
                                --singularity-args "--bind $mntpt:$mntpt" \
                                --use-singularity \
                                --snakefile $filename
                            ;;
                        yes)
                            nohup snakemake \
                                  --configfile config/${HOSTNAME}.yaml \
                                  --cores $threads \
                                  --directory ${repo} \
                                  --singularity-args "--bind $mntpt:$mntpt" \
                                  --use-singularity \
                                  --snakefile $filename
                            ;;
                    esac
                done                
                ;;
            force_final)
                source activate snakemake
                select nohup_option in no yes
                do
                    case $nohup_option in
                        no)                                          
                            snakemake \
                                --configfile config/${HOSTNAME}.yaml \
                                --cores $threads \
                                --directory ${repo} \
                                --force \
                                --singularity-args "--bind $mntpt:$mntpt" \
                                --use-singularity \
                                --snakefile $filename
                            ;;
                        yes)
                            nohup snakemake \
                                  --configfile config/${HOSTNAME}.yaml \
                                  --cores $threads \
                                  --directory ${repo} \
                                  --force \
                                  --singularity-args "--bind $mntpt:$mntpt" \
                                  --use-singularity \
                                  --snakefile $filename
                            ;;
                    esac
                done                
                ;;            
            force_all)
                source activate snakemake
                select nohup_option in no yes
                do
                    case $nohup_option in
                        no)                                          
                
                snakemake \
                    --configfile config/${HOSTNAME}.yaml \
                    --cores $threads \
                    --directory ${repo} \
                    -F \
                    --singularity-args "--bind $mntpt:$mntpt" \
                    --use-singularity \
                    --snakefile $filename
                                            ;;
                        yes)
                nohup snakemake \
                    --configfile config/${HOSTNAME}.yaml \
                    --cores $threads \
                    --directory ${repo} \
                    -F \
                    --singularity-args "--bind $mntpt:$mntpt" \
                    --use-singularity \
                    --snakefile $filename
                            ;;
                    esac
                done                
                ;;            
        esac
        break
    done
    break
done
