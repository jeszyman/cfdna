# [[file:~/repos/cfdna-wgs/cfdna-wgs.org::*Adapter-trim%20and%20QC%20reads%20with%20fastp][Adapter-trim and QC reads with fastp:2]]
#!/usr/bin/env bash
set -o errexit   # abort on nonzero exitstatus
set -o nounset   # abort on unbound variable
set -o pipefail  # don't hide errors within pipes

# Script variables

input_read1="${1}"
input_read2="${2}"
log_html="${3}"
log_json="${4}"
output_read1="${5}"
output_read2="${6}"
output_failed="${7}"
output_unpaired1="${8}"
output_unpaired2="${9}"
params_threads="${10}"

# Functions
main(){
    fastp_wrap $output_failed \
               $input_read1 \
               $input_read2 \
               $log_html \
               $log_json \
               $output_read1 \
               $output_read2 \
               $output_unpaired1 \
               $output_unpaired2 \
               $params_threads
}

fastp_wrap(){
    fastp --detect_adapter_for_pe \
          --failed_out $output_failed \
          --in1 $input_read1 \
          --in2 $input_read2 \
          --html $log_html \
          --json $log_json \
          --out1 $output_read1 \
          --out2 $output_read2 \
          --unpaired1 $output_unpaired1 \
          --unpaired2 $output_unpaired2 \
          --thread $params_threads
    }

# Run
main "$@"
# Adapter-trim and QC reads with fastp:2 ends here
