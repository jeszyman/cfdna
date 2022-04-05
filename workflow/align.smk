configfile: "./config/common.yaml"

rule all:
    input:
        config["inputs_dir"] + "/" + config["hg38_fasta"],
        config["inputs_dir"] + "/" + config["hg38_bwa_index_zip"],
        config["ref_dir"] + "/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.amb",
        config["ref_dir"] + "/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.ann",
        config["ref_dir"] + "/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwt",
        config["ref_dir"] + "/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.sa",

rule get_fasta:
    log: config["log_dir"] + "/get_fasta.log"
    output:
        hg38_fasta = config["inputs_dir"] + "/" + config["hg38_fasta"],
        hg38_bwa_index_zip = config["inputs_dir"] + "/" + config["hg38_bwa_index_zip"],
    shell:
        """
        wget {config[hg38_fasta_ftp]} -P {config[inputs_dir]}
        wget {config[hg38_bwa_index_ftp]} -P {config[inputs_dir]}
        tar -xzf {config[inputs_dir]}/{config[hg38_bwa_index_zip]} -C {config[ref_dir]}
        """
