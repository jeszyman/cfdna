# [[file:~/repos/cfdna-wgs/cfdna-wgs.org::*Make%20alignment%20index][Make alignment index:2]]
#!/usr/bin/env bash
in_fasta="${1}"
out_prefix="${2}"

bwa index -p $out_prefix $in_fasta
# Make alignment index:2 ends here
