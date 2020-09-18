// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process SUBREAD_FEATURECOUNTS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "quay.io/biocontainers/subread:2.0.1--hed695b0_0"
    //container "https://depot.galaxyproject.org/singularity/subread:2.0.1--hed695b0_0"

    conda (params.conda ? "bioconda::subread=2.0.1" : null)

    input:
    tuple val(meta), path(bams), path(annotation)
    val   options

    output:
    tuple val(meta), path("*featureCounts.txt")        , emit: counts
    tuple val(meta), path("*featureCounts.txt.summary"), emit: summary
    path "*.version.txt"                               , emit: version

    script:
    def software   = getSoftwareName(task.process)
    def ioptions   = initOptions(options)
    def prefix     = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"
    def paired_end = meta.single_end ? '' : '-p'

    def strandedness = 0
    if (meta.strandedness == 'forward') {
        strandedness = 1
    } else if (meta.strandedness == 'reverse') {
        strandedness = 2
    }
    """
    featureCounts \\
        $ioptions.args \\
        $paired_end \\
        -T $task.cpus \\
        -a $annotation \\
        -s $strandedness \\
        -o ${prefix}.featureCounts.txt \\
        ${bams.join(' ')}

    echo \$(featureCounts -v 2>&1) | sed -e "s/featureCounts v//g" > ${software}.version.txt
    """
}
