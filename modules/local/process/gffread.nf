// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process GFFREAD {
    tag "$gff"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda     (params.enable_conda ? "bioconda::gffread=0.12.1" : null)
    container "quay.io/biocontainers/gffread:0.12.1--h8b12597_0"

    input:
    path fasta
    path gff
    
    output:
    path "*.fa"         , emit: fasta
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    gffread $options.args -F -w transcripts.fa -g $fasta $gff
    echo \$(gffread --version 2>&1) > ${software}.version.txt
    """
}
