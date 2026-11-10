// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process CUSTOM_MULTIQCCUSTOMBIOTYPE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6b/6bb45b3a3fed2ebd55af1322a927ea11a751b8b62e542021dbf5be0d375abbb7/data'
        : 'community.wave.seqera.io/library/python:3.12.12--2411dfb0536aecdf'}"

    input:
    tuple val(meta), path(count)
    tuple val(meta2), path(header)

    output:
    tuple val(meta), path("*biotype_counts_mqc.tsv")      , emit: tsv
    tuple val(meta), path("*biotype_counts_rrna_mqc.tsv") , emit: rrna
    path "versions.yml"                                   , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'mqc_features_stat.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.biotype_counts_mqc.tsv
    touch ${prefix}.biotype_counts_rrna_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
