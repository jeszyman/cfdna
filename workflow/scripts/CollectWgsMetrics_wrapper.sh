input=$1
picard_jar=$2
genome_fasta=$3
output=$4

java -jar $picard_jar CollectWgsMetrics \
       INPUT=input \
       OUTPUT=output \
       READ_LENGTH=150 \
       REFERENCE_SEQUENCE=genome_fasta
