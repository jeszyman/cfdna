#########1#########2#########3#########4#########5#########6#########7#########8
#                                                                              #
#     Fragmentomic Analysis of Cell-free DNA Whole Genome Sequencing           #
#                                                                              #
#########1#########2#########3#########4#########5#########6#########7#########8

rule make_gc_map_bind:
    container: frag_container,
    input:
        gc5mb = config["gc5mb"],
        blklist = config["blklist"],
    log: logdir + "/make_gc_map_bind.log",
    output: refdir + "/keep_5mb.bed",
    params:
        script = "{frag_script_dir}/make_gc_map_bind.sh",
    shell:
        """
        {params.script} \
        {input.gc5mb} \
        {input.blklist} \
        {output} &> {log}
        """

# Make a bed file from filtered bam
rule filt_bam_to_frag_bed:
    benchmark: benchdir + "/{library}_filt_bam_to_frag_bed.benchmark.txt",
    container: frag_container,
    input: frag_bams + "/{library}_filt.bam",
    log: logdir + "/{library}_filt_bam_to_frag_bed.log",
    output: frag_frag_beds + "/{library}_filt.bed",
    params:
        fasta = genome_fasta,
        script = "{frag_script_dir}/filt_bam_to_frag_bed.sh",
        threads = frag_threads,
    shell:
        """
        {params.script} \
	{input} \
        {params.fasta} \
        {params.threads} \
        {output}
        """

# Make GC distributions
rule gc_distro:
    benchmark: benchdir + "/{library}_frag_gc_distro.benchmark.txt",
    container: frag_container,
    input: frag_frag_beds + "/{library}_filt.bed",
    log: logdir + "/{library}_frag_gc_distro.log",
    output: frag_frag_gc_distros + "/{library}_gc_distro.csv",
    params:
        script = "{frag_script_dir}/gc_distro.R",
    shell:
        """
        Rscript {params.script} \
        {input} \
        {output} \
        > {log} 2>&1
        """

# Make healthy GC distributions summary file
rule healthy_gc:
    benchmark: benchdir + "/frag_healthy_gc.benchmark.txt",
    container: frag_container,
    input: expand(frag_frag_gc_distros + "/{library}_gc_distro.csv", library = FRAG_HEALTHY_LIBRARIES),
    log: logdir + "/frag_healthy_gc.log",
    output: frag_frag_gc_distros + "/healthy_med.rds",
    params:
        distro_dir = frag_frag_gc_distros,
        script = "{frag_script_dir}/healthy_gc.R",
    shell:
        """
        Rscript {params.script} \
        {params.distro_dir} \
        "{input}" \
        {output} > {log} 2>&1
        """

# Sample fragments by healthy GC proportions
rule frag_gc_sample:
    benchmark: benchdir + "/{library}_frag_gc_sample.benchmark.txt",
    container: frag_container,
    input:
        frag_bed = frag_frag_beds + "/{library}_filt.bed",
        healthy_med = frag_frag_gc_distros + "/healthy_med.rds",
    log: logdir + "/{library}_frag_gc_sample.log",
    output: frag_frag_beds + "/{library}_sampled_frag.bed",
    params:
        script = "{frag_script_dir}/gc_sample.R",
    shell:
        """
        Rscript {params.script} \
        {input.healthy_med} \
        {input.frag_bed} \
        {output} > {log} 2>&1
        """

# Sum fragments in short and long length groups

rule frag_sum:
    benchmark: benchdir + "/{library}_frag_sum.benchmark.txt",
    container: frag_container,
    input: frag_frag_beds + "/{library}_sampled_frag.bed",
    log: logdir + "/{library}_frag_frag_window_sum.log",
    output:
        short = frag_frag_beds + "/{library}_norm_short.bed",
        long =  frag_frag_beds + "/{library}_norm_long.bed",
    params:
        script = "{frag_script_dir}/frag_window_sum.sh",
        threads = frag_threads,
    shell:
        """
        {params.script} \
        {input} \
        {output.short} {output.long} &> {log}
        """

# Count short and long fragments intersecting kept genomic windows

rule frag_window_count:
    benchmark: benchdir + "/{library}_frag_frag_window_int.benchmark.txt",
    container: frag_container,
    input:
        short = frag_frag_beds + "/{library}_norm_short.bed",
        long = frag_frag_beds + "/{library}_norm_long.bed",
        matbed = refdir + "/keep_5mb.bed",
    log: logdir + "/{library}_frag_frag_window_int.log",
    output:
        short = frag_frag_counts + "/{library}_cnt_short.tmp",
        long = frag_frag_counts + "/{library}_cnt_long.tmp",
    params:
        script = "{frag_script_dir}/frag_window_int.sh",
        threads = threads,
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
rule frag_count_merge:
    benchmark: benchdir + "/frag_count_merge.benchmark.txt",
    container: frag_container,
    input: expand(frag_frag_counts + "/{library}_cnt_{length}.tmp",  library = FRAG_LIBS, length = ["short","long"]),
    log: logdir + "/frag_count_merge.log",
    output:  frag_frag + "/frag_counts.tsv",
    params:
        counts_dir = frag_frag + "/counts",
        script = "{frag_script_dir}/count_merge.sh",
        threads = frag_threads,
    shell:
        """
        {params.script} \
        {params.counts_dir} \
        {output} &> {log}
        """

rule unit_cent_sd:
    benchmark: benchdir + "/unit_cent_sd.benchmark.txt",
    container: frag_container,
    input: frag_frag + "/frag_counts.tsv",
    log: logdir + "/unit_cent_sd.log",
    output: frag_frag + "/ratios.tsv",
    params:
        script = "{frag_script_dir}/make_ratios.R",
    shell:
        """
        Rscript {params.script} \
        {input} {output} > {log} 2>&1
        """
