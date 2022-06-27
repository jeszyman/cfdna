#!/usr/bin/env bash
input=$1
threads=$2
blacklist=$3
output=$4

bamPEFragmentSize --bamfiles $input \
                  --numberOfProcessors $threads \
                  --blackListFileName $blacklist \
                  --maxFragmentLength 500 \
                  --outRawFragmentLengths $output
