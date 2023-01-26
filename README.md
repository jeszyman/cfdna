
# Introduction

This repository hosts a snakemake workflow for basic processing of whole-genome sequencing reads from cell-free DNA.

![img](resources/int_test.png)


# Organization

Master branch of the repository contains most recent developments while stable versions are saved as terminal branches (*e.g.* stable.1.0.0).

Directory `workflow` contains two types of workflows- process-focused snakefiles (reads.smk, cna.smk, frag.smk) suitable for integration into another snakemake pipeline using the :include command, and the \_int\_test snakefile with examples of such integration using the repository test data.


# Use

-   All software needed for the pipeline is present within the associated docker container (see `docker` and <https://hub.docker.com/repository/docker/jeszyman/cfdna_wgs/general>).
-   See the example configuration yaml `config/int_test.yaml` and wrapper workflow `workflow/int_test` for necessary run conditions.


# Changelog

-   <span class="timestamp-wrapper"><span class="timestamp">[2023-01-26 Thu] </span></span> - Version 9.0.0: Removed -f 3 flag for perfectly matched pairs in samtools filtering as the black from BWA removes some fragments at a set max length. Added framework for benchmark analysis. Added conditional execution of downsampling. Removed (temporarily) final wig and ichor commands of CNA as these don't currently run correctly without full genome alignment, so can't be validated on test data. Added local documentation of cfdna-wgs dockerfile.
-   <span class="timestamp-wrapper"><span class="timestamp">[2023-01-21 Sat] </span></span> - Version 8.0.0: Corrected rule filt\_bam\_to\_frag\_bed to fix mates of inputs, which seems to prevent errors in the bamtobed call. Frag\_window\_count now uses windows of consistent 5 Mb size, which are generated from rule make\_gc\_map\_bind. Added a merged fragment counts file and zero-centered unit SD counts.
-   <span class="timestamp-wrapper"><span class="timestamp">[2022-12-07 Wed] </span></span> - Version 7.0.0: Added copy number alteration and DELFI fragmentomics.
-   <span class="timestamp-wrapper"><span class="timestamp">[2022-10-17 Mon] </span></span> - Version 6.0.0: Using fastp for read trimming (replaces trimmomatic). Simplified naming schema. Removed downsampling (will reinstate in later version).
-   <span class="timestamp-wrapper"><span class="timestamp">[2022-09-08 Thu] </span></span> - Version 5.3.0: some minor name changes
-   <span class="timestamp-wrapper"><span class="timestamp">[2022-08-19 Fri] </span></span> - Version 5.2.0 validated: Adds bamCoverage and plotCoverage from deeptools. Benchmarks BWA.
-   <span class="timestamp-wrapper"><span class="timestamp">[2022-08-09 Tue] </span></span> - Version 5.1.0 validated: Added cfdna wgs-specific container for each rule, referenced to config
-   <span class="timestamp-wrapper"><span class="timestamp">[2022-08-05 Fri] </span></span> - Version 5.0.0 validated: Added a symlink rule based on python dictionary. Added repo-specific output naming, added checks for sequence type and file readability to input tsv.
-   <span class="timestamp-wrapper"><span class="timestamp">[2022-06-27 Mon] </span></span> - Version 4 validated. Further expanded read\_qc.tsv table. Removed bam post-processing step and added a more expansive bam filtering step. Updated downsampling to work off filtered alignments.
-   <span class="timestamp-wrapper"><span class="timestamp">[2022-06-26 Sun] </span></span> - Version 3.2 validated. Expanded the qc aggregate table and added some comments.
-   <span class="timestamp-wrapper"><span class="timestamp">[2022-06-24 Fri] </span></span> - Validate version 3.1 which includes genome index build as a snakefile rule.
-   <span class="timestamp-wrapper"><span class="timestamp">[2022-06-24 Fri] </span></span> - Validated version 3 with read number checkpoint for down-sampling.
-   <span class="timestamp-wrapper"><span class="timestamp">[2022-05-31 Tue] </span></span> - Conforms to current biotools best practices.
-   <span class="timestamp-wrapper"><span class="timestamp">[2022-04-29 Fri] </span></span> - Moved multiqc to integration testing as inputs are dependent on final sample labels. Integration testing works per this commit.
