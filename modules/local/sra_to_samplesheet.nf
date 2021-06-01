// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

params.options     = [:]
params.results_dir = ''

process SRA_TO_SAMPLESHEET {
    tag "$meta.id"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    memory 100.MB

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*csv"), emit: csv

    exec:
    //  Remove custom keys needed to download the data
    def meta_map = meta.clone()
    meta_map.remove("id")
    meta_map.remove("fastq_1")
    meta_map.remove("fastq_2")
    meta_map.remove("md5_1")
    meta_map.remove("md5_2")
    meta_map.remove("single_end")

    // Add required fields for the pipeline to the beginning of the map
    pipeline_map = [
        sample      : "${meta.id.split('_')[0..-2].join('_')}",
        fastq_1     : "${params.outdir}/${params.results_dir}/${fastq[0]}",
        fastq_2     : meta.single_end ? '' : "${params.outdir}/${params.results_dir}/${fastq[1]}",
        strandedness: 'unstranded'
    ]
    pipeline_map << meta_map

    // Create a samplesheet
    csv  = pipeline_map.keySet().collect{ '"' + it + '"'}.join(",") + '\n'
    csv += pipeline_map.values().collect{ '"' + it + '"'}.join(",")

    // Write to file
    def file = task.workDir.resolve("${meta.id}.samplesheet.csv")
    file.text = csv
}
