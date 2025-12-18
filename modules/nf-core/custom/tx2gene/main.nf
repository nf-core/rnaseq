process CUSTOM_TX2GENE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10.4' :
        'biocontainers/python:3.10.4' }"

    input:
    tuple val(meta), path(gtf)
    tuple val(meta2), path ("quants/*")
    val quant_type
    val id
    val extra

    output:
    tuple val(meta), path("*tx2gene.tsv"), emit: tx2gene
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'tx2gene.py'

    stub:
    """
    touch ${meta.id}.tx2gene.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
