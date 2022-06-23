#!/usr/bin/env bash

# Remove output directories
outputs_dirs=("fastq"
              "bam"
              "logs"
              "qc")

for dir in "${outputs_dirs[@]}"; do
    if [ -d test/${dir} ]; then \rm -rf test/${dir}; fi
done
