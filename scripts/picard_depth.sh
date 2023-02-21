# [[file:~/repos/cfdna-wgs/cfdna-wgs.org::*Sequencing%20depth%20metrics%20via%20Picard][Sequencing depth metrics via Picard:2]]
#!/usr/bin/env bash
input=$1
picard_jar=$2
genome_fasta=$3
output=$4

java -jar $picard_jar CollectWgsMetrics \
       INPUT=$input \
       OUTPUT=$output \
       READ_LENGTH=150 \
       REFERENCE_SEQUENCE=$genome_fasta
# Sequencing depth metrics via Picard:2 ends here
