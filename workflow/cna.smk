#########1#########2#########3#########4#########5#########6#########7#########8
#                                                                              #
#                    Copy-number Alteration Analysis of                        #
#                  Cell-free DNA Whole Genome Sequencing                       #
#                                                                              #
#                                                                              #
#########1#########2#########3#########4#########5#########6#########7#########8

# Filter fragments by length
rule cna_frag_filt:
    benchmark: benchdir + "/{library}_{frag_distro}_cfdna_wgs_frag_filt.benchmark.txt",
    container: cfdna_wgs_container,
    input: cfdna_wgs_cna_in_bams + "/{library}.bam",
    log: logdir + "/{library}_{frag_distro}_cfdna_wgs_frag_filt.log",
    output:
        nohead = temp(cfdna_wgs_cna_frag_bams) + "/{library}_frag{frag_distro}.nohead",
        onlyhead = temp(cfdna_wgs_cna_frag_bams) + "/{library}_frag{frag_distro}.only",
        final = cfdna_wgs_cna_frag_bams + "/{library}_frag{frag_distro}.bam",
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
rule bam_to_wig:
    benchmark: benchdir + "/{library}_{frag_distro}_cfdna_wgs_bam_to_wig.benchmark.txt",
    container: cfdna_wgs_container,
    input: cfdna_wgs_cna_frag_bams + "/{library}_frag{frag_distro}.bam",
    log: logdir + "/{library}_{frag_distro}_cfdna_wgs_bam_to_wig.log",
    output: cfdna_wgs_cna_wigs + "/{library}_frag{frag_distro}.wig",
    params:
        chrs = chrs,
        script = cfdna_wgs_scriptdir + "/bam_to_wig.sh",
        threads = cfdna_wgs_threads,
    shell:
        """
        /opt/hmmcopy_utils/bin/readCounter \
        --chromosome "{params.chrs}" \
        --quality 20 \
        --window 1000000 \
        {input} > {output}
        """

# Run ichorCNA without a panel of normals
rule ichor_nopon:
    input:
        wig = cfdna_wgs_cna_wigs + "/{library}_frag{frag_distro}.wig",
    output:
        cfdna_wgs_cna_ichor_nopon + "/{library}_frag{frag_distro}.cna.seg",
    params:
        script = cfdna_wgs_scriptdir + "/MOD_runIchorCNA.R",
        out_dir = cfdna_wgs_cna_ichor_nopon,
    container:
        cfdna_wgs_container,
    shell:
        """
        Rscript {params.script} \
         --id {wildcards.library}_frag{wildcards.frag_distro} \
         --WIG {input.wig} \
         --gcWig /opt/ichorCNA/inst/extdata/gc_hg38_1000kb.wig \
         --mapWig /opt/ichorCNA/inst/extdata/map_hg38_1000kb.wig \
         --centromere /opt/ichorCNA/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt \
         --normal "c(0.95, 0.99, 0.995, 0.999)" \
         --ploidy "c(2)" \
         --maxCN 3 \
         --estimateScPrevalence FALSE \
         --scStates "c()" \
         --outDir {params.out_dir}
        """
