#!/usr/bin/env bash
outputs_dirs=("fastq"
              "bam"
              "logs"
              "qc")

for dir in "${outputs_dirs[@]}"; do
    if [ -d test/${dir} ]; then \rm -rf test/${dir}; fi
done

../basecamp/src/smk_dry_run.sh config/int_test.yaml workflow/int_test.smk &&
../basecamp/src/smk_draw.sh config/int_test.yaml workflow/int_test.smk resources/int_test.pdf &&
../basecamp/src/smk_forced_run.sh config/int_test.yaml workflow/int_test.smk &&
echo "Integration testing passed, do you want to erase results files?" &&
select yn in "Yes" "No"; do
    case $yn in
        Yes ) for dir in "${outputs_dirs[@]}"; do \rm -rf test/$dir; done; break;;
        No ) exit;;
    esac
done
