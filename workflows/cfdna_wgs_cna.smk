
rule cfdna_wgs_fastp:
    conda:
        "../config/cfdna-wgs-conda-env.yaml"
    input:
        r1 = f"{data_dir}/cfdna-wgs/fastqs/{{library_id}}.raw_R1.fastq.gz",
        r2 = f"{data_dir}/cfdna-wgs/fastqs/{{library_id}}.raw_R2.fastq.gz"
    log:
        html = f"{data_dir}/logs/{{library_id}}_cfdna_wgs_fastp.html",
        json = f"{data_dir}/logs/{{library_id}}_cfdna_wgs_fastp.json",
    output:
        failed = f"{data_dir}/cfdna-wgs/fastqs/{{library_id}}.failed.fastq.gz",
        r1 = f"{data_dir}/cfdna-wgs/fastqs/{{library_id}}.processed_R1.fastq.gz",
        r2 = f"{data_dir}/cfdna-wgs/fastqs/{{library_id}}.processed_R2.fastq.gz",
        up1 = f"{data_dir}/cfdna-wgs/fastqs/{{library_id}}.unpaired_R1.fastq.gz",
        up2 = f"{data_dir}/cfdna-wgs/fastqs/{{library_id}}.unpaired_R2.fastq.gz",
    shell:
        """
        fastp --detect_adapter_for_pe \
        --failed_out {output.failed} \
        --in1 {input.r1} --in2 {input.r2} \
        --html {log.html} --json {log.json} \
        --out1 {output.r1} --out2 {output.r2} \
        --unpaired1 {output.up1} --unpaired2 {output.up2} \
        --thread 8
        """
rule cfdna_wgs_fastqc:
    conda:
        "../config/cfdna-wgs-conda-env.yaml"
    input:
        f"{data_dir}/cfdna-wgs/fastqs/{{library_id}}.{{processing}}_{{read}}.fastq.gz",
    log:
        f"{data_dir}/logs/{{library_id}}.{{processing}}_{{read}}_cfdna_wgs_fastqc.log",
    output:
        f"{data_dir}/cfdna-wgs/qc/{{library_id}}.{{processing}}_{{read}}_fastqc.html",
        f"{data_dir}/cfdna-wgs/qc/{{library_id}}.{{processing}}_{{read}}_fastqc.zip",
    params:
        outdir = f"{data_dir}/cfdna-wgs/qc",
        threads = 2,
    resources:
        concurrency = 25,
    shell:
        """
        fastqc \
        --outdir {params.outdir} \
        --quiet \
        --svg \
        --threads {params.threads} \
        {input} &> {log}
        """
rule cfdna_wgs_bwa_index:
    conda:
        "../config/cfdna-wgs-conda-env.yaml"
    input:
        lambda wildcards: f"{data_dir}/inputs/{config['cfdna_wgs_ref_assemblies'][wildcards.name]['input']}",
    output:
        fa = f"{data_dir}/ref/bwa/{{name}}/{{name}}.fa",
        amb = f"{data_dir}/ref/bwa/{{name}}/{{name}}.fa.amb",
        ann = f"{data_dir}/ref/bwa/{{name}}/{{name}}.fa.ann",
        bwt = f"{data_dir}/ref/bwa/{{name}}/{{name}}.fa.bwt",
        pac = f"{data_dir}/ref/bwa/{{name}}/{{name}}.fa.pac",
        sa  = f"{data_dir}/ref/bwa/{{name}}/{{name}}.fa.sa",
    params:
        fasta_target = lambda wildcards: f"{data_dir}/ref/bwa/{wildcards.name}/{wildcards.name}.fa",
    log:
        f"{data_dir}/logs/{{name}}_bwa_index.log"
    shell:
        """
        mkdir -p $(dirname {params.fasta_target})
        cp {input} {params.fasta_target}
        bwa index -p {params.fasta_target[:-3]} {params.fasta_target} > {log} 2>&1
        """
rule cfdna_wgs_bwa_mem:
    conda:
        "../config/cfdna-wgs/conda-env.yaml"
    input:
        r1 = f"{data_dir}/analysis/cfdna-wgs/fastqs/{{library_id}}_processed_R1.fastq.gz",
        r2 = f"{data_dir}/analysis/cfdna-wgs/fastqs/{{library_id}}_processed_R2.fastq.gz",
        ref = f"{data_dir}/ref/bwa/{{ref_name}}/{{ref_name}}.fa",
        sa_check = f"{data_dir}/ref/bwa/{{ref_name}}/{{ref_name}}.fa.sa",
    output:
        bam = f"{data_dir}/analysis/cfdna-wgs/bams/{{library_id}}.{{ref_name}}.bwa.coorsort.bam",
    shell:
        """
        bwa mem -M -t 24 \
        {input.ref} {input.r1} {input.r2} \
        | samtools view -@ 4 -Sb - -o - \
        | samtools sort -@ 4 - -o {output.bam}
        samtools index -@ 4 {output.bam}
