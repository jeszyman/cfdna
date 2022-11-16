# Cell-free DNA whole genome sequencing analysis of copy number alteration

# Filter fragments by length
rule cfdna_wgs_frag_filt:
    benchmark: logdir + "/{library}_{milreads}_{frag_distro}_cfdna_wgs_frag_filt.benchmark.txt",
    container: cfdna_wgs_container,
    input: analysis + "/cfdna_wgs_bams/{library}_ds{milreads}.bam",
    log: logdir + "/{library}_{milreads}_{frag_distro}_cfdna_wgs_frag_filt.log",
    output:
        nohead = temp(analysis + "/cfdna_wgs_frag/{library}_ds{milreads}_frag{frag_distro}.nohead"),
        onlyhead = temp(analysis + "/cfdna_wgs_frag/{library}_ds{milreads}_frag{frag_distro}.only"),
        final = analysis + "/cfdna_wgs_frag/{library}_ds{milreads}_frag{frag_distro}.bam",
    params:
        script = cfdna_wgs_scriptdir + "/frag_filt.sh",
        threads = cfdna_wgs_threads,
    shell:
        """
        frag_min=$(echo {wildcards.frag_distro} | sed -e "s/_.*$//g")
        frag_max=$(echo {wildcards.frag_distro} | sed -e "s/^.*_//g")
        {params.script} \
        {input} \
        {output.nohead} \
        $frag_min \
        $frag_max \
        {config[threads]} \
        {output.onlyhead} \
        {output.final}
        """

# Use readCounter to create windowed wig from bam file
rule cfdna_wgs_bam_to_wig:
    benchmark: logdir + "/{library}_{milreads}_{frag_distro}_cfdna_wgs_bam_to_wig.benchmark.txt",
    container: cfdna_wgs_container,
    input: analysis + "/cfdna_wgs_frag/{library}_ds{milreads}_frag{frag_distro}.bam",
    log: logdir + "/{library}_{milreads}_{frag_distro}_cfdna_wgs_bam_to_wig.log",
    output: analysis + "/cfdna_wgs_frag/{library}_ds{milreads}_frag{frag_distro}.wig",
    params:
        chrs = chrs,
        script = cfdna_wgs_scriptdir + "/bam_to_wig.sh",
        threads = cfdna_wgs_threads,
    shell:
        """
        /opt/hmmcopy_utils/bin/readCounter --window 1000000 --quality 20 \
        --chromosome {params.chrs} \
        {input} > {output} &> {log}
        """
