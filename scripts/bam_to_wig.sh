#!/usr/bin/env bash
input=$1
output=$2

        /opt/hmmcopy_utils/bin/readCounter --window 1000000 --quality 20 \
        --chromosome {params.chrs} \
        {input} > {output}
