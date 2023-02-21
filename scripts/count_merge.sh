# [[file:~/repos/cfdna-wgs/cfdna-wgs.org::*Merge%20counts%20across%20length%20and%20library][Merge counts across length and library:2]]
# For unit testing
#counts_dir="/home/jeszyman/mpnst/analysis/cfdna-wgs/frag/counts"
#out_tsv="/home/jeszyman/mpnst/analysis/cfdna-wgs/frag/frag_counts.tsv"

# Define variables
counts_dir="${1}"
out_tsv="${2}"

# Remove the existing aggregate file if present
if [ -f $out_tsv ]; then rm $out_tsv; fi
#touch $out_tsv

# Make aggregate file
for file in ${counts_dir}/*;
do
    # Add file name to each line
    awk '{{print FILENAME (NF?"\t":"") $0}}' $file |
        # Modify file name to library id
        sed 's/^.*lib/lib/g' |
        sed 's/_.*_/\t/g' |
        # Cleanup "tmp"
        sed 's/.tmp//g' |
        # Send to output
        sed 's/\.bed//g' >> $out_tsv
done

# Add a header
sed -i  '1 i\library	len_class	chr	start	end	gc	count' $out_tsv
# Merge counts across length and library:2 ends here
