process CUSTOM_CATADDITIONALFASTA {
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(fasta), path(gtf)
    tuple val(meta2), path(add_fasta)
    val(biotype)

    output:
    tuple val(meta), path("out/${prefix}.fasta"), emit: fasta
    tuple val(meta), path("out/${prefix}.gtf")  , emit: gtf
    path "versions.yml"                         , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    template 'fasta2gtf.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir out
    touch out/${prefix}.fasta
    touch out/${prefix}.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
