#!/usr/bin/env bash
filelist=$1
out_dir=$2

Rscript /opt/ichorCNA/scripts/createPanelOfNormals.R --filelist $filelist \
        --chrs "paste0('chr', c(1:22, \"X\"))" \
        --chrNormalize "c(1:22, \"X\")" \
        --gcWig /opt/ichorCNA/inst/extdata/gc_hg38_1000kb.wig \
        --mapWig /opt/ichorCNA/inst/extdata/map_hg38_1000kb.wig \
        --centromere /opt/ichorCNA/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt  \
        --outfile "${out_dir}/pon"
