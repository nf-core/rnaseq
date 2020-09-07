// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process QUALIMAP_RNASEQ {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "quay.io/biocontainers/qualimap:2.2.2d--1"
    //container "https://depot.galaxyproject.org/singularity/qualimap:2.2.2d--1"

    conda (params.conda ? "bioconda::qualimap=2.2.2d" : null)

    input:
    tuple val(meta), path(bam)
    path  gtf
    val   options

    output:
    tuple val(meta), path("${prefix}"), emit: results
    path  "*.version.txt"             , emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    prefix       = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"

    // Figure out strandedness from pipeline parameters
    def unstranded       = params.unstranded
    def forward_stranded = params.forward_stranded
    def reverse_stranded = params.reverse_stranded
    if (params.pico) {
        unstranded       = false
        forward_stranded = true
        reverse_stranded = false
    }
    def strandedness = 'non-strand-specific'
    if (forward_stranded) {
        strandedness = 'strand-specific-forward'
    } else if (reverse_stranded) {
        strandedness = 'strand-specific-reverse'
    }
    def paired_end = meta.single_end ? '' : '-pe'
    def memory     = task.memory.toGiga() + "G"
    """
    unset DISPLAY
    mkdir tmp
    export _JAVA_OPTIONS=-Djava.io.tmpdir=./tmp
    qualimap \\
        --java-mem-size=$memory \\
        rnaseq \\
        $ioptions.args \\
        -bam $bam \\
        -gtf $gtf \\
        -p $strandedness \\
        $paired_end \\
        -outdir $prefix

    echo \$(qualimap 2>&1) | sed 's/^.*QualiMap v.//; s/Built.*\$//' > ${software}.version.txt
    """
}
