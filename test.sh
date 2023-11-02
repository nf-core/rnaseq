#!/usr/bin/env bash

nextflow run . \
    --outdir output \
    --pseudo_aligner kallisto \
    -resume -profile conda,test
