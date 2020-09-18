// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process PICARD_MARKDUPLICATES {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "quay.io/biocontainers/picard:2.23.2--0"
    //container "https://depot.galaxyproject.org/singularity/picard:2.23.2--0"

    conda (params.conda ? "bioconda::picard=2.23.2" : null)

    input:
    tuple val(meta), path(bam)
    val   options

    output:
    tuple val(meta), path("*.bam")        , emit: bam
    tuple val(meta), path("*.metrics.txt"), emit: metrics
    path  "*.version.txt"                 , emit: version

    script:
    def software  = getSoftwareName(task.process)
    def ioptions  = initOptions(options)
    def prefix    = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"
    def avail_mem = 3
    if (!task.memory) {
        log.info '[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    picard \\
        -Xmx${avail_mem}g \\
        MarkDuplicates \\
        $ioptions.args \\
        INPUT=$bam \\
        OUTPUT=${prefix}.bam \\
        METRICS_FILE=${prefix}.MarkDuplicates.metrics.txt

    echo \$(picard MarkDuplicates --version 2>&1) | awk -F' ' '{print \$NF}' > ${software}.version.txt
    """
}
