# Aggregate QC files using MultiQC
rule frag_multiqc:
    input:
        lambda wildcards: expand(f"{qc_dir}/{{library}}_cfdna_wgs_fastp.json",
                                 library = lib_map[wildcards.lib_set]['libs']),
        lambda wildcards: expand(f"{qc_dir}/{{library}}_{{processing}}_{{read}}_fastqc.zip",
                                 library = lib_map[wildcards.lib_set]['libs'],
                                 processing = ['raw','proc',],
                                 read = ['R1','R2']),
        lambda wildcards: expand(f"{qc_dir}/{{library}}_{{build}}_{{processing}}_samstats.txt",
                                 library = lib_map[wildcards.lib_set]['libs'],
                                 build = lib_map[wildcards.lib_set]['build'],
                                 processing = ['raw','dedup','filt']),
        lambda wildcards: expand(f"{qc_dir}/{{library}}_{{build}}_{{processing}}_flagstat.txt",
                                 library = lib_map[wildcards.lib_set]['libs'],
                                 build = lib_map[wildcards.lib_set]['build'],
                                 processing = ['raw','dedup','filt']),
        lambda wildcards: expand(f"{qc_dir}/{{library}}_{{build}}_picard_depth.txt",
                                 library = lib_map[wildcards.lib_set]['libs'],
                                 build = lib_map[wildcards.lib_set]['build']),
        f"{qc_dir}/deeptools_{{lib_set}}_lengths.txt",
        f"{qc_dir}/{{lib_set}}_frag_coverage.tsv",
    log: f"{log_dir}/{{lib_set}}_cfdna_wgs_multiqc.log",
    output:
        f"{qc_dir}/{{lib_set}}_cfdna_wgs_multiqc.html",
        f"{qc_dir}/{{lib_set}}_cfdna_wgs_multiqc_fastqc.txt",
        f"{qc_dir}/{{lib_set}}_cfdna_wgs_multiqc_data/multiqc_samtools_stats.txt",
        f"{qc_dir}/{{lib_set}}_cfdna_wgs_multiqc_data/multiqc_picard_wgsmetrics.txt",
        f"{qc_dir}/{{lib_set}}_cfdna_wgs_multiqc_data/multiqc_samtools_flagstat.txt",
    params:
        out_dir = f"{qc_dir}",
        out_name = "frag_multiqc",
        script = f"{cfdna_script_dir}/multiqc.sh",
        threads = threads,
    shell:
        """
        {params.script} \
        "{input}" \
        {params.out_name} \
        {params.out_dir} &> {log}
        """
# Make a tab-separated aggregate QC table
rule make_cfdna_wgs_qc_tsv:
    input:
        fq = f"{qc_dir}/{{lib_set}}_cfdna_wgs_multiqc_data/multiqc_fastqc.txt",
        mqsam = f"{qc_dir}/{{lib_set}}_cfdna_wgs_multiqc_data/multiqc_samtools_stats.txt",
        mqflag = f"{qc_dir}/{{lib_set}}_cfdna_multiqc_data/multiqc_samtools_flagstat.txt",
        picard = f"{qc_dir}/{{lib_set}}_multiqc_data/multiqc_picard_wgsmetrics.txt",
        deeptools_frag = f"{qc_dir}/{{lib_set}}_deeptools_cfdna_wgs_lengths.txt",
        deeptools_cov = f"{qc_dir}/{{lib_set}}_cfdna_wgs_coverage.tsv",
    log: f"{log_dir}/{{lib_set}}_cfdna_wgs_make_qc_tsv.log",
    output:
        readqc = f"{qc_dir}/{{lib_set}}_cfdna_wgs_read_qc.tsv",
        fraglen = f"{qc_dir}/{{lib_set}}_cfdna_wgs_len.tsv",
    params:
        script = f"{cfdna_script_dir}/make_qc_tsv.R",
    shell:
        """
        Rscript {params.script} \
        {input.fq} \
        {input.mqsam} \
        {input.mqflag} \
        {input.picard} \
        {input.deeptools_frag} \
        {input.deeptools_cov} \
        {output.readqc} \
        {output.fraglen} >& {log}
        """
