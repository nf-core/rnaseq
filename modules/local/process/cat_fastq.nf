// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]
def options    = initOptions(params.options)

/*
 * Concatenate FastQ files
 */
process CAT_FASTQ {
    tag "$meta.id"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'merged_fastq', publish_id:meta.id) }

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }

    // conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    // } else {
    //     container "biocontainers/biocontainers:v1.2.0_cv1"
    // }
    
    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("*.merged.fastq.gz"), emit: reads

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    readList     = reads.collect{it.toString()}
    if (!meta.single_end) {
        if (readList.size > 2) {
            def read1 = []
            def read2 = []
            readList.eachWithIndex{ v, ix -> ( ix & 1 ? read2 : read1 ) << v }
            """
            cat ${read1.sort().join(' ')} > ${prefix}_1.merged.fastq.gz
            cat ${read2.sort().join(' ')} > ${prefix}_2.merged.fastq.gz
            """
        } else {
            """
            ln -s ${reads[0]} ${prefix}_1.merged.fastq.gz
            ln -s ${reads[1]} ${prefix}_2.merged.fastq.gz
            """
        }
    } else {
        if (readList.size > 1) {
            """
            cat ${readList.sort().join(' ')} > ${prefix}.merged.fastq.gz
            """
        } else {
            """
            ln -s $reads ${prefix}.merged.fastq.gz
            """
        }
    }
}
