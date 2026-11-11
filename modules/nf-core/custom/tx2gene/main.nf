// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process CUSTOM_TX2GENE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a3/a3afdbbeb3c6f212ac595eca51a566b047a7bab2c50783955d3aac0d55fa683b/data'
        : 'community.wave.seqera.io/library/python:3.10.4--b5b45c494792cf95'}"

    input:
    tuple val(meta), path(gtf)
    tuple val(meta2), path ("quants/*")
    val quant_type
    val id
    val extra

    output:
    tuple val(meta), path("*tx2gene.tsv"), emit: tx2gene
    path "versions.yml"                  , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'tx2gene.py'

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.tx2gene.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
