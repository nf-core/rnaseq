// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process EAUTILS_GTF2BED {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a9/a90490336208590dfdaaf102bb64f3c341aa0ced1f95d9247020ce88ae017db7/data'
        : 'community.wave.seqera.io/library/perl:5.26.2--2f879f377ee6109b'}"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("${prefix}.bed") , emit: bed
    path "versions.yml"                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    args   = task.ext.args ?: ''

    """
    echo $args
    """

    template 'gtf2bed.pl'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl -e 'print substr(\$^V, 1)')
    END_VERSIONS
    """
}
