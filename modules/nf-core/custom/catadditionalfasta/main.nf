process CUSTOM_CATADDITIONALFASTA {
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(fasta), path(gtf)
    tuple val(meta2), path(add_fasta)
    val  biotype

    output:
    tuple val(meta), path("*/*.fasta") , emit: fasta
    tuple val(meta), path("*/*.gtf")   , emit: gtf
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'fasta2gtf.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir out
    touch out/genome_transcriptome.fasta
    touch out/genome_transcriptome.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | grep -v "Python ")
    END_VERSIONS
    """
}
