// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process GFFREAD {
    tag "$gff"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container "quay.io/biocontainers/gffread:0.11.7--h8b12597_0"
    //container https://depot.galaxyproject.org/singularity/gffread:0.11.7--h8b12597_0

    conda (params.conda ? "bioconda::gffread=0.11.7" : null)

    input:
    path fasta
    path gff
    val  options

    output:
    path "*.fa"         , emit: fasta
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    """
    gffread $ioptions.args -F -w transcripts.fa -g $fasta $gff
    echo \$(gffread --version 2>&1) > ${software}.version.txt
    """
}
