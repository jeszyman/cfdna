#!/usr/bin/env bash

fastqc --outdir $2 \
       --quiet \
       --threads $3 $1
