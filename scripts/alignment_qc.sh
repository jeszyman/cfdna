# [[file:~/repos/cfdna-wgs/cfdna-wgs.org::*Alignment%20QC][Alignment QC:2]]
#!/usr/bin/env bash
set -o errexit   # abort on nonzero exitstatus
set -o nounset   # abort on unbound variable
set -o pipefail  # don't hide errors within pipes

# Script variables
input="${1}"
log_flagstat="${2}"
log_samstat="${3}"
output_flagstat="${4}"
output_samstat="${5}"
threads="${6}"

# Functions
main(){
    flagstat $input $output_flagstat $log_flagstat $threads
    samstats $input $output_samstat $log_samstat $threads
}

flagstat(){
    local input="${1}"
    local output="${2}"
    local log="${3}"
    local threads="${4}"
    #
    samtools flagstat -@ $threads $input > $output 2>$log
}

samstats(){
    local input="${1}"
    local output="${2}"
    local log="${3}"
    local threads="${4}"
    #
    samtools stats -@ $threads $input > $output 2>$log
}

# Run
main "$@"
# Alignment QC:2 ends here
