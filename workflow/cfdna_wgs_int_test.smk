container: config["container"]

IDS, = glob_wildcards(config["fq_dir"] + "/{id}_R1.fastq.gz")
MILREADS = config["MILREADS"]
	   
rule all:
    input:
        expand(config["processed_fq_dir"] + "/{read_id}_proc_{read}.fastq.gz", read_id = IDS, read = ["R1","R2"]),
        expand(config["unpr_fq_dir"] + "/{read_id}_unpr_R1.fastq.gz", read_id = IDS, read = ["R1","R2"]),
        expand(config["bam_dir"] + "/{read_id}.sam", read_id = IDS),
        expand(config["qc_dir"] + "/{read_id}_{read}_fastqc.html", read_id = IDS, read = ["R1","R2"]),
        expand(config["qc_dir"] + "/{read_id}_proc_{read}_fastqc.html", read_id = IDS, read = ["R1","R2"]),
        expand(config["bam_dir"] + "/{read_id}_dedup.bam", read_id = IDS),
        expand(config["bam_dir"] + "/{read_id}_dedup.bam.bai", read_id = IDS),
        expand(config["qc_dir"] + "/{read_id}_{bam_step}_samstats.txt", read_id = IDS, bam_step= ["dedup","raw"]),
        expand(config["qc_dir"] + "/{read_id}_{bam_step}_flagstat.txt", read_id = IDS, bam_step =["dedup","raw"]),
        expand(config["bam_dir"] + "/{read_id}_ds{milreads}.bam", read_id = IDS, milreads = MILREADS),
        config["qc_dir"] + "/all_qc.html",

include: "read_preprocess.smk"

rule multiqc:
    input:
        expand(config["qc_dir"] + "/{read_id}_{read}_fastqc.html", read_id = IDS, read = ["R1","R2"]),
        expand(config["qc_dir"] + "/{read_id}_proc_{read}_fastqc.html", read_id = IDS, read = ["R1","R2"]),
        expand(config["qc_dir"] + "/{read_id}_{bam_step}_samstats.txt", read_id = IDS, bam_step= ["dedup","raw"]),
        expand(config["qc_dir"] + "/{read_id}_{bam_step}_flagstat.txt", read_id = IDS, bam_step =["dedup","raw"]),
    params:
        out_dir = config["qc_dir"]
    output:
        config["qc_dir"] + "/all_qc.html"
    shell:
        """
        multiqc {params.out_dir} \
        --force \
        --outdir {params.out_dir} \
        --filename all_qc 
        """
