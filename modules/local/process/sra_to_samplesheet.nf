// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

params.options     = [:]
params.results_dir = ''

/*
 * Stage FastQ files downloaded by SRA and auto-create a samplesheet for the pipeline
 */
process SRA_TO_SAMPLESHEET {
    tag "$meta.id"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

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
        group       : "${meta.id.split('_')[0..-2].join('_')}",
        replicate   : 1,
        fastq_1     : [ "${params.outdir}", "${params.results_dir}", "${fastq[0]}" ].join(File.separator),
        fastq_2     : meta.single_end ? '' : [ "${params.outdir}", "${params.results_dir}", "${fastq[1]}" ].join(File.separator),
        strandedness: 'unstranded'
    ]
    pipeline_map << meta_map

    // Write to file
    File file = new File(["${task.workDir}", "${meta.id}.samplesheet.csv"].join(File.separator))
    file.write pipeline_map.keySet().collect{ '"' + it + '"'}.join(",") + '\n'
    file.append(pipeline_map.values().collect{ '"' + it + '"'}.join(",")) + '\n'
}