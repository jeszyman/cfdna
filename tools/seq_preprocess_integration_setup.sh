#!/bin/echo Run:.

# For documentation, not intended to be executable

if [ -d test ]; then \rm -rf test; fi
mkdir -p test/fastq
zcat /mnt/ris/aadel/mpnst/inputs/MPNST/19_2_082_R1.fastq.gz | head -n 100000 > "test/fastq/mpnst1_R1.fastq"
zcat /mnt/ris/aadel/mpnst/inputs/MPNST/19_2_082_R2.fastq.gz | head -n 100000 > "test/fastq/mpnst1_R2.fastq"
zcat /mnt/ris/aadel/mpnst/inputs/MPNST/25_2_072_R1.fastq.gz | head -n 100000 > "test/fastq/mpnst2_R1.fastq"
zcat /mnt/ris/aadel/mpnst/inputs/MPNST/25_2_072_R2.fastq.gz | head -n 100000 > "test/fastq/mpnst2_R2.fastq"
zcat /mnt/ris/aadel/mpnst/inputs/PN/37_JS0050CD112717_R1.fastq.gz | head -n 100000 > "test/fastq/plex1_R1.fastq"
zcat /mnt/ris/aadel/mpnst/inputs/PN/37_JS0050CD112717_R2.fastq.gz | head -n 100000 > "test/fastq/plex1_R2.fastq"
zcat /mnt/ris/aadel/mpnst/inputs/PN/30_JS0044CD112818_R1.fastq.gz | head -n 100000 > "test/fastq/plex2_R1.fastq"
zcat /mnt/ris/aadel/mpnst/inputs/PN/30_JS0044CD112818_R2.fastq.gz | head -n 100000 > "test/fastq/plex2_R2.fastq"
for file in "test/fastq/*.fastq"; do gzip $file; done

mkdir -p "test/inputs"
wget --directory-prefix="test/inputs/" https://raw.githubusercontent.com/usadellab/Trimmomatic/main/adapters/TruSeq3-PE.fa
wget --directory-prefix="test/inputs/" https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

cp resources/samples.tsv test/inputs/

mkdir -p test/ref
zcat "test/inputs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz" | grep -A 2000 chr8 > test/inputs/chr8.fa
\rm test/inputs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

singularity shell ~/sing_containers/biotools.sif
bwa index -p test/ref/chr8 test/inputs/chr8.fa
exit
