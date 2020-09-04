// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

/*
 * Concatenate FastQ files
 */
process CAT_FASTQ {
    tag "$meta.id"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.conda ? "${baseDir}/environment.yml" : null)

    input:
    tuple val(meta), path(reads)
    val   options

    output:
    tuple val(meta), path("*.merged.fastq.gz"), emit: reads

    script:
    def ioptions = initOptions(options)
    def prefix   = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"
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
