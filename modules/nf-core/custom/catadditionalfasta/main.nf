// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process CUSTOM_CATADDITIONALFASTA {
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/dd/dd0e85ebfdd196217d90b852b82b98f50e1f222f43aa55119358e558e7a8a71f/data'
        : 'community.wave.seqera.io/library/python:3.12.2--b229f007be7a84c0'}"

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
