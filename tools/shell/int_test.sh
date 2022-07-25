#!/usr/bin/env bash
shopt -s extglob
cd test
\rm -rf !(inputs)
cd ../

smk_dry_run.sh config/int_test.yaml workflow/int_test.smk \
    && smk_draw.sh config/int_test.yaml workflow/int_test.smk resources/int_test.pdf \
    && smk_forced_run.sh config/int_test.yaml workflow/int_test.smk \
    && echo "Integration testing passed, do you want to erase results files?" \
    && select yn in "Yes" "No"; do
           case $yn in
               Yes )
                   shopt -s extglob
                   cd test
                   \rm -rf !(inputs)
                   cd ../; break;;
               No ) exit;;
           esac
       done
