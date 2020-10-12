// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process STRINGTIE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "quay.io/biocontainers/stringtie:2.1.4--h7e0af3c_0"
    //container "https://depot.galaxyproject.org/singularity/stringtie:2.1.4--h7e0af3c_0"

    conda (params.conda ? "bioconda::stringtie=2.1.4" : null)

    input:
    tuple val(meta), path(bam)
    path  gtf
    val   options

    output:
    tuple val(meta), path("*.coverage.gtf")   , emit: coverage_gtf
    tuple val(meta), path("*.transcripts.gtf"), emit: transcript_gtf
    tuple val(meta), path("*.txt")            , emit: abundance
    tuple val(meta), path("*.ballgown")       , emit: ballgown
    path  "*.version.txt"                     , emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    def prefix   = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"

    def strandedness = ''
    if (meta.strandedness == 'forward') {
        strandedness = '--fr'
    } else if (meta.strandedness == 'reverse') {
        strandedness = '--rf'
    }
    """
    stringtie \\
        $bam \\
        $strandedness \\
        -G $gtf \\
        -o ${prefix}.transcripts.gtf \\
        -A ${prefix}.gene_abundance.txt \\
        -C ${prefix}.coverage.gtf \\
        -b ${prefix}.ballgown \\
        $ioptions.args

    stringtie --version > ${software}.version.txt
    """
}
