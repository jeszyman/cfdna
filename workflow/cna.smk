#########1#########2#########3#########4#########5#########6#########7#########8
#                                                                              #
#                    Copy-number Alteration Analysis of                        #
#                  Cell-free DNA Whole Genome Sequencing                       #
#                                                                              #
#                                                                              #
#########1#########2#########3#########4#########5#########6#########7#########8

# Use readCounter to create windowed wig from bam file
rule bam_to_wig:
    benchmark: benchdir + "/{library}_ds{downsample}_{frag_distro}_cfdna_wgs_bam_to_wig.benchmark.txt",
    container: cfdna_wgs_container,
    input: cfdna_wgs_bams + "/{library}_ds{downsample}_frag{frag_distro}.bam",
    log: logdir + "/{library}_ds{downsample}_{frag_distro}_cfdna_wgs_bam_to_wig.log",
    output: cfdna_wgs_wigs + "/{library}_ds{downsample}_frag{frag_distro}.wig",
    params:
        chrs = chrs,
        outdir = cfdna_wgs_wigs,
        script = cfdna_wgs_scriptdir + "/bam_to_wig.sh",
        threads = cfdna_wgs_threads,
    shell:
        """
        mkdir -p {params.outdir}
        /opt/hmmcopy_utils/bin/readCounter \
        --chromosome "{params.chrs}" \
        --quality 20 \
        --window 1000000 \
        {input} > {output}
        """

# Run ichorCNA without a panel of normals
rule ichor_nopon:
    input: cfdna_wgs_wigs + "/{library}_ds{downsample}_frag{frag_distro}.wig",
    output: cfdna_wgs_ichor_nopon + "/{library}_ds{downsample}_frag{frag_distro}.cna.seg",
    params:
        script = cfdna_wgs_scriptdir + "/MOD_runIchorCNA.R",
        out_dir = cfdna_wgs_ichor_nopon,
    container:
        cfdna_wgs_container,
    shell:
        """
        Rscript {params.script} \
         --id {wildcards.library}_frag{wildcards.frag_distro} \
         --WIG {input} \
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
