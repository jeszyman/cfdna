# Cell-free DNA whole genome sequencing fragmentomics

# Make a bed file from filtered bam
rule cfdna_wgs_filt_bam_to_frag_bed:
    benchmark: logdir + "/{library}_cfdna_wgs_filt_bam_to_frag_bed.benchmark.txt",
    container: cfdna_wgs_container,
    input: cfdna_wgs_frag_input_bams + "/{library}.bam",
    log: logdir + "/{library}_cfdna_wgs_filt_bam_to_frag_bed.log",
    output: cfdna_wgs_frag_beds + "/{library}_filt.bed",
    params:
        fasta = genome_fasta,
        script = cfdna_wgs_scriptdir + "/filt_bam_to_frag_bed.sh",
        threads = cfdna_wgs_threads,
    resources: mem_mb=5000
    shell:
        """
        {params.script} \
	{input} \
        {params.fasta} \
        {params.threads} \
        {output}
        """

# Make GC distributions
rule cfdna_wgs_gc_distro:
    benchmark: logdir + "/{library}_cfdna_wgs_gc_distro.benchmark.txt",
    container: cfdna_wgs_container,
    input: cfdna_wgs_frag_beds + "/{library}_filt.bed",
    log: logdir + "/{library}_cfdna_wgs_gc_distro.log",
    output: cfdna_wgs_frag_gc_distros + "/{library}_gc_distro.csv",
    params:
        script = cfdna_wgs_scriptdir + "/gc_distro.R",
    shell:
        """
        Rscript {params.script} \
        {input} \
        {output} \
        > {log} 2>&1
        """

# Make healthy GC distributions summary file
rule cfdna_wgs_healthy_gc:
    benchmark: logdir + "/cfdna_wgs_healthy_gc.benchmark.txt",
    container: cfdna_wgs_container,
    input: expand(cfdna_wgs_frag_gc_distros + "/{library}_gc_distro.csv", library = CFDNA_WGS_HEALTHY_LIBRARIES),
    log: logdir + "/cfdna_wgs_healthy_gc.log",
    output: cfdna_wgs_frag_gc_distros + "/healthy_med.rds",
    params:
        distro_dir = cfdna_wgs_frag_gc_distros,
        script = cfdna_wgs_scriptdir + "/healthy_gc.R",
    shell:
        """
        Rscript {params.script} \
        {params.distro_dir} \
        "{input}" \
        {output} \
        > {log} 2>&1
        """

# Sample fragments by healthy GC proportions
rule cfdna_wgs_gc_sample:
    benchmark: logdir + "/{library}_cfdna_wgs_gc_sample.benchmark.txt",
    container: cfdna_wgs_container,
    input:
        frag_bed = cfdna_wgs_frag_beds + "/{library}_filt.bed",
        healthy_med = cfdna_wgs_frag_gc_distros + "/healthy_med.rds",
    log: logdir + "/{library}_cfdna_wgs_gc_sample.log",
    output: cfdna_wgs_frag_beds + "/{library}_sampled_frag.bed",
    params:
        script = cfdna_wgs_scriptdir + "/gc_sample.R",
    shell:
        """
        Rscript {params.script} \
        {input.healthy_med} \
        {input.frag_bed} \
        {output} \
        > {log} 2>&1
        """

# Sum fragments in genomic windows by length
rule cfdna_wgs_frag_window_sum:
    benchmark: logdir + "/{library}_cfdna_wgs_frag_window_sum.benchmark.txt",
    container: cfdna_wgs_container,
    input: cfdna_wgs_frag_beds + "/{library}_sampled_frag.bed",
    log: logdir + "/{library}_cfdna_wgs_frag_window_sum.log",
    output:
        short = cfdna_wgs_frag_beds + "/{library}_norm_short.bed",
        long =  cfdna_wgs_frag_beds + "/{library}_norm_long.bed",
    params:
        script = cfdna_wgs_scriptdir + "/frag_window_sum.sh",
        threads = cfdna_wgs_threads,
    shell:
        """
        {params.script} \
        {input} \
        {output.short} {output.long} &> {log}
        """

# Count fragments intersecting windows
rule cfdna_wgs_frag_window_int:
    benchmark: logdir + "/{library}_cfdna_wgs_frag_window_int.benchmark.txt",
    container: cfdna_wgs_container,
    input:
        short = cfdna_wgs_frag_beds + "/{library}_norm_short.bed",
        long = cfdna_wgs_frag_beds + "/{library}_norm_long.bed",
        matbed = keep_bed,
    log: logdir + "/{library}_cfdna_wgs_frag_window_int.log",
    output:
        short = cfdna_wgs_frag_counts + "/{library}_cnt_short.tmp",
        long = cfdna_wgs_frag_counts + "/{library}_cnt_long.tmp",
    params:
        script = cfdna_wgs_scriptdir + "/frag_window_int.sh",
        threads = cfdna_wgs_threads,
    shell:
        """
        {params.script} \
        {input.short} \
        {input.matbed} \
        {output.short}
        {params.script} \
        {input.long} \
        {input.matbed} \
        {output.long}
        """

# Merge short and long fragment counts by genomic poistion for all libraries
rule cfdna_wgs_count_merge:
    benchmark: logdir + "/cfdna_wgs_count_merge.benchmark.txt",
    container: cfdna_wgs_container,
    input: expand(cfdna_wgs_frag_counts + "/{library}_cnt_{length}.tmp",  library = CNA_WGS_LIBRARIES, length = ["short","long"]),
    log: logdir + "/cfdna_wgs_count_merge.log",
    output:  cfdna_wgs_frag + "/frag_counts.tsv",
    params:
        script = cfdna_wgs_scriptdir + "/count_merge.sh",
        threads = cfdna_wgs_threads,
    shell:
        """
        array=({input})
        {params.script} \
        {output} "$array" &> {log}
        """
