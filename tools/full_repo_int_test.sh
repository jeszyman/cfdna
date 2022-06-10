# Remove results directories 
results_dirs=()
readarray -t results_dirs <<< $(ls -d1 test/* | grep -v -e '^test/f' -e ref -e inputs)
printf '%s\n' "${results_dirs[@]}"


for dir in "${results_dirs[@]}"; 
do
\rm -rf $dir
done

# Run test 
basecamp/src/smk_dry_run.sh config/int_repo_test.yaml workflow/cfdna_wgs_int_test.smk
basecamp/src/smk_forced_run.sh config/int_repo_test.yaml workflow/cfdna_wgs_int_test.smk
