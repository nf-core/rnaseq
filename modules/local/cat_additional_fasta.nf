// Import generic module functions
include { saveFiles; getProcessName } from './functions'

params.options = [:]

process CAT_ADDITIONAL_FASTA {
    tag "$add_fasta"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'genome', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }

    input:
    path fasta
    path gtf
    path add_fasta
    val  biotype

    output:
    path "${name}.fasta", emit: fasta
    path "${name}.gtf"  , emit: gtf

    script:
    def genome_name  = params.genome ? params.genome : fasta.getBaseName()
    def biotype_name = biotype ? "-b $biotype" : ''
    def add_name     = add_fasta.getBaseName()
    name             = "${genome_name}_${add_name}"
    """
    fasta2gtf.py -o ${add_fasta.baseName}.gtf $biotype_name $add_fasta
    cat $fasta $add_fasta > ${name}.fasta
    cat $gtf ${add_fasta.baseName}.gtf > ${name}.gtf
    """
}
