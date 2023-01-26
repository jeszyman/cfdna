#########1#########2#########3#########4#########5#########6#########7#########8
#                                                                              #
#                    Basic Read and Alignment Processing of                    #
#                    Cell-free DNA Whole Genome Sequencing                     #
#                                                                              #
#########1#########2#########3#########4#########5#########6#########7#########8

# Make alignment index
#  Note: Upon first run, this rule will touch an empty file with the same path
#        as the index prefix. Thereafter, you can avoid repeat indexing when the
#        rule "sees" this empty file. For repo intergration testing with an
#        external reference, indexing can likewise be avoided with this empty
#        file at the external index location.

rule cfdna_wgs_index:
    benchmark: benchdir + "/cfdna_wgs_index.benchmark.txt",
    container: cfdna_wgs_container,
    input: genome_fasta,
    log: logdir + "/cfdna_wgs_index.log",
    output: done = touch(genome_ref),
    params:
        out_prefix = genome_ref,
        script = cfdna_wgs_scriptdir + "/index.sh",
        threads = cfdna_wgs_threads,
    shell:
        """
        bwa index -p {params.out_prefix} {input} &> {log}
        """

# Adapter-trim and QC reads with fastp
rule cfdna_wgs_fastp:
    benchmark: benchdir + "/{library}_cfdna_wgs_fastp.benchmark.txt",
    container: cfdna_wgs_container,
    input:
        read1 = cfdna_wgs_fastqs + "/{library}_raw_R1.fastq.gz",
        read2 = cfdna_wgs_fastqs + "/{library}_raw_R2.fastq.gz",
    log:
        cmd = logdir + "/{library}_cfdna_wgs_fastp.log",
        html = logdir + "/{library}_cfdna_wgs_fastp.html",
        json = logdir + "/{library}_cfdna_wgs_fastp.json",
    output:
        read1 = cfdna_wgs_fastqs + "/{library}_processed_R1.fastq.gz",
        read2 = cfdna_wgs_fastqs + "/{library}_processed_R2.fastq.gz",
        failed = cfdna_wgs_fastqs + "/{library}_failed_fastp.fastq.gz",
        unpaired1 = cfdna_wgs_fastqs + "/{library}_unpaired_R1.fastq.gz",
        unpaired2 = cfdna_wgs_fastqs + "/{library}_unpaired_R2.fastq.gz",
    params:
        script = cfdna_wgs_scriptdir + "/fastp.sh",
        threads = cfdna_wgs_threads,
    resources:
        mem_mb = 500,
    shell:
        """
        {params.script} \
        {input.read1} \
        {input.read2} \
        {log.html} \
        {log.json} \
        {output.read1} \
        {output.read2} \
        {output.failed} \
        {output.unpaired1} \
        {output.unpaired2} \
        {params.threads} &> {log.cmd}
        """

# Align reads with BWA
rule cfdna_wgs_align:
    benchmark: benchdir + "/{library}_cfdna_wgs_align.benchmark.txt",
    container: cfdna_wgs_container,
    input:
        ref = genome_ref,
        read1 = cfdna_wgs_fastqs + "/{library}_processed_R1.fastq.gz",
        read2 = cfdna_wgs_fastqs + "/{library}_processed_R2.fastq.gz",
    log: logdir + "/{library}_cfdna_wgs_align.log",
    output:
        sort = cfdna_wgs_bams + "/{library}_raw.bam",
        index = cfdna_wgs_bams + "/{library}_raw.bam.bai",
    params:
        script = cfdna_wgs_scriptdir + "/align.sh",
        threads = 4,
    resources:
        mem_mb = 500,
    shell:
        """
        {params.script} \
        {input.ref} \
        {input.read1} \
        {input.read2} \
        {params.threads} \
        {output.sort} &> {log}
        """

# Remove PCR duplicates from aligned reads
rule cfdna_wgs_dedup:
    benchmark: benchdir + "/{library}_cfdna_wgs_dedup.benchmark.txt",
    container: cfdna_wgs_container,
    input: cfdna_wgs_bams + "/{library}_raw.bam",
    log: logdir + "/{library}_cfdna_wgs_dedup.log",
    output: cfdna_wgs_bams + "/{library}_dedup.bam",
    params:
        script = cfdna_wgs_scriptdir + "/dedup.sh",
        threads = cfdna_wgs_threads,
    shell:
        """
        {params.script} \
        {input} \
        {output} \
        {params.threads} &> {log}
        """

# Filter de-duplicated alignments.
# Remove unmapped, not primary, and duplicate reads. Additional location filter by config bedfile variable.

checkpoint cfdna_wgs_filter_alignment:
    benchmark: benchdir + "/{library}_cfdna_wgs_filter_alignment.benchmark.txt",
    container: cfdna_wgs_container,
    input: cfdna_wgs_bams + "/{library}_dedup.bam",
    log: logdir + "/{library}_cfdna_wgs_filter_alignment.log",
    output: cfdna_wgs_bams + "/{library}_filt.bam",
    params:
        script = cfdna_wgs_scriptdir + "/filter_alignment.sh",
        threads = cfdna_wgs_threads,
    shell:
        """
        {params.script} \
        {input} \
        {params.threads} \
        {output} &> {log}
        """

# Get read quality by FASTQC
rule cfdna_wgs_fastqc:
    benchmark: benchdir+ "/{library}_{processing}_{read}_cfdna_wgs_fastqc.benchmark.txt",
    container: cfdna_wgs_container,
    input: cfdna_wgs_fastqs + "/{library}_{processing}_{read}.fastq.gz",
    log: logdir + "/{library}_{processing}_{read}_cfdna_wgs_fastqc.log",
    output:
        qcdir + "/{library}_{processing}_{read}_fastqc.html",
        qcdir + "/{library}_{processing}_{read}_fastqc.zip",
    params:
        outdir = qcdir,
        script = cfdna_wgs_scriptdir + "/fastqc.sh",
        threads = cfdna_wgs_threads,
    shell:
        """
        {params.script} \
        {input} \
        {params.outdir} \
        {params.threads} &> {log}
        """

# Get alignment QC using samtools
rule cfdna_wgs_alignment_qc:
    container: cfdna_wgs_container,
    input: cfdna_wgs_bams + "/{library}_{processing}.bam",
    log:
        flagstat = logdir + "/{library}_{processing}_flagstat_cfdna_wgs_alignment_qc.log",
        samstat = logdir + "/{library}_{processing}_samstats_cfdna_wgs_alignment_qc.log",
    output:
        flagstat = qcdir + "/{library}_{processing}_flagstat.txt",
        samstat = qcdir + "/{library}_{processing}_samstats.txt",
    params:
        script = cfdna_wgs_scriptdir + "/alignment_qc.sh",
        threads = cfdna_wgs_threads,
    shell:
        """
        {params.script} \
        {input} \
        {log.flagstat} \
        {log.samstat} \
        {output.flagstat} \
        {output.samstat} \
        {params.threads}
        """

# Sequencing depth metrics via Picard
rule cfdna_wgs_picard_depth:
    benchmark: benchdir + "/{library}_cfdna_wgs_picard_depth.benchmark.txt",
    container: cfdna_wgs_container,
    input: cfdna_wgs_bams + "/{library}_filt.bam",
    log: logdir + "/{library}_cfdna_wgs_picard_depth.log",
    output: qcdir + "/{library}_picard_depth.txt",
    params:
        script = cfdna_wgs_scriptdir + "/picard_depth.sh",
        threads = cfdna_wgs_threads,
    shell:
        """
        {params.script} \
        {input} \
        {config[picard_jar]} \
        {config[genome_fasta]} \
        {output}
        """

# Get fragment sizes using deepTools
rule cfdna_wgs_bampefragsize:
    benchmark: benchdir + "/cfdna_wgs_bampefragsize.benchmark.txt",
    container: cfdna_wgs_container,
    input: expand(cfdna_wgs_bams + "/{library}_filt.bam", library = CFDNA_WGS_LIBRARIES),
    log: logdir + "/cfdna_wgs_bampefragsize.log",
    output:
        raw = qcdir + "/deeptools_frag_lengths.txt",
        hist = qcdir + "/deeptools_frag_lengths.png",
    params:
        blacklist = config["blklist"],
        script = cfdna_wgs_scriptdir + "/bampefragsize.sh",
        threads = cfdna_wgs_threads,
    shell:
        """
        {params.script} \
        "{input}" \
        {log} \
        {output.hist} \
        {output.raw} \
        {params.blacklist} \
        {params.threads}
        """

# Make deeptools bamCoverage bedfile
rule cfdna_wgs_bamcoverage:
    benchmark: benchdir + "/{library}_cfdna_wgs_bamcoverage.benchmark.txt",
    container: cfdna_wgs_container,
    input: cfdna_wgs_bams + "/{library}_filt.bam",
    log: logdir + "/{library}_cfdna_wgs_bamcoverage.log",
    output: qcdir + "/{library}_bamcoverage.bg",
    params:
        bin = "10000",
        blacklist = config["blklist"],
        script = cfdna_wgs_scriptdir + "/bamcoverage.sh",
        threads = cfdna_wgs_threads,
    shell:
        """
        {params.script} \
        {input} \
        {output} \
        {params.bin} \
        {params.blacklist} \
        {params.threads} &> {log}
        """

# Make deepTools plotCoverage coverage maps for all filtered bams
rule cfdna_wgs_plotcoverage:
    benchmark: benchdir + "/cfdna_wgs_plotcoverage.benchmark.txt",
    container: cfdna_wgs_container,
    input: expand(cfdna_wgs_bams + "/{library}_filt.bam", library = CFDNA_WGS_LIBRARIES),
    log: logdir + "/cfdna_wgs_plotcoverage.log",
    output:
        raw = qcdir + "/cfdna_wgs_coverage.tsv",
        plot = qcdir + "/cfdna_wgs_coverage.pdf",
    params:
        blacklist = config["blklist"],
        script = cfdna_wgs_scriptdir + "/plotcoverage.sh",
        threads = cfdna_wgs_threads,
    shell:
        """
        {params.script} \
        "{input}" \
        {params.blacklist} \
        {params.threads} \
        {output.raw} \
        {output.plot} &> {log}
        """

# Aggregate QC files using MultiQC
rule cfdna_wgs_multiqc:
    benchmark: benchdir + "/cfdna_wgs_multiqc.benchmark.txt",
    container: cfdna_wgs_container,
    input:
        expand(logdir + "/{library}_cfdna_wgs_fastp.json", library = CFDNA_WGS_LIBRARIES),
        expand(qcdir + "/{library}_{processing}_{read}_fastqc.zip", library = CFDNA_WGS_LIBRARIES, processing = ["raw", "processed", "unpaired"], read = ["R1","R2"]),
        expand(qcdir + "/{library}_{processing}_samstats.txt", library = CFDNA_WGS_LIBRARIES, processing = ["raw","filt"]),
        expand(qcdir + "/{library}_{processing}_flagstat.txt", library = CFDNA_WGS_LIBRARIES, processing = ["raw","filt"]),
        expand(qcdir + "/{library}_picard_depth.txt", library = CFDNA_WGS_LIBRARIES),
        qcdir + "/deeptools_frag_lengths.txt",
        qcdir + "/cfdna_wgs_coverage.tsv",
    log: logdir + "/cfdna_wgs_multiqc.log",
    output:
        qcdir + "/cfdna_wgs_multiqc.html",
        qcdir + "/cfdna_wgs_multiqc_data/multiqc_fastqc.txt",
        qcdir + "/cfdna_wgs_multiqc_data/multiqc_samtools_stats.txt",
        qcdir + "/cfdna_wgs_multiqc_data/multiqc_picard_wgsmetrics.txt",
        qcdir + "/cfdna_wgs_multiqc_data/multiqc_samtools_flagstat.txt",
    params:
        out_dir = qcdir,
        out_name = "cfdna_wgs_multiqc",
        script = cfdna_wgs_scriptdir + "/multiqc.sh",
        threads = cfdna_wgs_threads,
    shell:
        """
        {params.script} \
        "{input}" \
        {params.out_name} \
        {params.out_dir} &> {log}
        """

# Make a tab-separated aggregate QC table
checkpoint cfdna_wgs_make_qc_tsv:
    benchmark: benchdir + "/cfdna_wgs_make_qc_tsv.benchmark.txt",
    container: cfdna_wgs_container,
    input:
        fq = qcdir + "/cfdna_wgs_multiqc_data/multiqc_fastqc.txt",
        mqsam = qcdir + "/cfdna_wgs_multiqc_data/multiqc_samtools_stats.txt",
        mqflag = qcdir + "/cfdna_wgs_multiqc_data/multiqc_samtools_flagstat.txt",
        picard = qcdir + "/cfdna_wgs_multiqc_data/multiqc_picard_wgsmetrics.txt",
        deeptools_frag = qcdir + "/deeptools_frag_lengths.txt",
        deeptools_cov = qcdir + "/cfdna_wgs_coverage.tsv",
    log: logdir + "/cfdna_wgs_make_qc_tsv.log",
    output:
        readqc = qcdir + "/cfdna_wgs_read_qc.tsv",
        fraglen = qcdir + "/cfdna_wgs_frag_len.tsv",
    params:
        script = cfdna_wgs_scriptdir + "/make_qc_tsv.R",
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

rule downsample_bams:
    container: cfdna_wgs_container,
    input: cfdna_wgs_bams + "/{library}_filt.bam",
    output: touch(logdir + "/{library}_{downsample}_downsample.done"),
    params:
        out_dir = cfdna_wgs_bams,
        script = cfdna_wgs_scriptdir + "/downsample_bams.sh",
        suffix = "_filt.bam",
        threads = cfdna_wgs_threads,
    shell:
        """
        {params.script} \
        {input} \
        {wildcards.downsample} \
        {params.out_dir} \
        {params.suffix} \
        {params.threads}
        """

# If downsample occured, then write filename into this per-library log, else leave the log file blank
rule log_dowsample:
    input: logdir + "/{library}_{downsample}_downsample.done",
    output: logdir + "/{library}_{downsample}_made",
    params:
        bamdir = cfdna_wgs_bams,
    shell:
        """
        dspath={params.bamdir}/{wildcards.library}_ds{wildcards.downsample}.bam
        if [ -f $dspath ]; then echo "$dspath"  > {output}; else touch {output}; fi
        """

# Use the downsampled bam logs to make a single text file of conditionally executed final targets.
# Specifically in this example, log text lines are in the form
# cfdna_wgs_bams + "/{library}_ds{downsample}_frag90_150.bam" to setup conditional execution of fragment filtering ONLY on downsampled bams
# Note alternative delimiter "~" to sed allows cfdna_wgs_wigs as param

checkpoint ds_cond_target_list:
    input: expand(logdir + "/{library}_{downsample}_made", library = CFDNA_WGS_LIBRARIES, downsample = DOWNSAMPLE),
    output: logdir + "/ds_final_targets",
    params:
        outdir = cfdna_wgs_bams,
        frag_distro=config["frag_distro"]
    shell:
        """
        if [ -f {output} ]; then rm {output}; fi
        cat {input} > {output}
        sed -i 's~^.*lib~{params.outdir}/lib~g' {output}
        sed -i 's/.bam$/_frag{params.frag_distro}.bam/g' {output}
        """

# Function jsut pulls the final target names out of ds_final_targets
def get_ds_targets(wildcards):
    with open(checkpoints.ds_cond_target_list.get(**wildcards).output[0], "r") as f:
      non_empty_files = [l.strip() for l in f.readlines()]
    return non_empty_files

# This rule allows execution of rules which will generate the conditional targets in ds_cond_target_list
rule make_ds_targets:
    input:
        get_ds_targets
    output: logdir + "/aggregate_output"
    run:
        with open(output[0], "w") as f:
            f.write("\n".join(input))

rule frag_filt:
    container: cfdna_wgs_container,
    input:
        main = cfdna_wgs_bams + "/{library}_ds{downsample}.bam",
        check = logdir + "/{library}_{downsample}_made",
    output:
        nohead = temp(cfdna_wgs_bams + "/{library}_ds{downsample}_frag{frag_distro}.nohead"),
        onlyhead = temp(cfdna_wgs_bams + "/{library}_ds{downsample}_frag{frag_distro}.only"),
        final = cfdna_wgs_bams + "/{library}_ds{downsample}_frag{frag_distro}.bam",
    params:
        script = cfdna_wgs_scriptdir + "/frag_filt.sh",
        threads = cfdna_wgs_threads,
    shell:
        """
        frag_min=$(echo {wildcards.frag_distro} | sed -e "s/_.*$//g")
        frag_max=$(echo {wildcards.frag_distro} | sed -e "s/^.*_//g")
        {params.script} \
        {input.main} \
        {output.nohead} \
        $frag_min \
        $frag_max \
        {config[threads]} \
        {output.onlyhead} \
        {output.final}
        """
