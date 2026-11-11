// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process CUSTOM_GTFFILTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/13/1340bf9bd68b2760272dfa0f7d62b3d6703505fd6899fd98407f2225e3917213/data'
        : 'community.wave.seqera.io/library/python:3.14.5--de7566bbce1bfa41'}"

    input:
    tuple val(meta), path(gtf)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("${prefix}.${suffix}"), emit: gtf
    path "versions.yml"                         , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    suffix = "gtf" + (gtf.extension == 'gz' ? '.gz' : '')
    args   = task.ext.args ?: ''

    """
    echo $args
    """

    template 'gtffilter.py'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    suffix = "gtf" + (gtf.extension == 'gz' ? '.gz' : '')
    """
    touch ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | cut -d ' ' -f 2)
    END_VERSIONS
    """
}
