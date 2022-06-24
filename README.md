This repository has a snakemake workflow for basic processing of whole-genome sequencing reads from cell-free DNA.

![img](resources/int_test.png)

Master branch of the repository contains most recent developments. Stable versions are saved as terminal branches (*e.g.* stable1.0.0).

Files labeled int\_test will run integration testing of all rules on a small dataset in test/inputs. See config/int\_test.yaml for necessary run conditions.


# Changlog

-   <span class="timestamp-wrapper"><span class="timestamp">[2022-06-24 Fri] </span></span> - Validated version 3 with read number checkpoint for down-sampling.
-   <span class="timestamp-wrapper"><span class="timestamp">[2022-05-31 Tue] </span></span> - Conforms to current biotools best practices.
-   <span class="timestamp-wrapper"><span class="timestamp">[2022-04-29 Fri] </span></span> - Moved multiqc to integration testing as inputs are dependent on final sample labels. Integration testing works per this commit.
