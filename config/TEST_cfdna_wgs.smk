ref_assemblies:
  ncbi_hg38:
    url:
      "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
    name: ncbi_hg38
    input: GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.gz

cfdna-cna-conda-env: "~/repos/cfdna/config/cfdna-cna-conda-env.yaml"
cfdna-cna-dir: "/mnt/data/projects/cfdna/cfdna-cna"
cfdna-scripts-dir: "~/repos/cfdna/scripts"
cfdna-repo-dir: "~/repos/cfdna"
ichor_repo: "~/repos/ichorCNA-patched"
local-tmp-dir: "/mnt/data/projects/cfdna/tmp"
log-dir: "/mnt/data/projects/cfdna/logs"

input-tsv: "~/repos/cfdna/data/TEST_cfdna_cna.tsv"
input-samples:
  - kh_01

samples-tsv: "~/repos/cfdna/data/TEST_cfdna_cna.tsv"

input-bams-dir: "/mnt/data/projects/cfdna/cfdna-cna/input_bams"

envs:
  cfdna-cna: "~/repos/cfdna/config/cfdna-cna-conda-env.yaml"

project-data-dir: "/mnt/data/projects/cfdna"

cfdna-repo: "~/repos/cfdna"
data-dir: "/mnt/data/projects/cfdna"
cfdna-wgs-dir: "/mnt/data/projects/cfdna/cfdna-wgs"
mosdepth-quant-levels: "1,5,10,20,30"
