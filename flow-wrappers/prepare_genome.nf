#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome.nf'

workflow {
    PREPARE_GENOME ( [

    ], "", false )
}