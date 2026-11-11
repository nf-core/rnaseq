// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process BOWTIE2_BUILD {
    tag "$fasta"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/70/7093a62c894ad3337cabda5efadce78a3ecdb5ba34cf8646743c6d43d8c4b77d/data'
:         'community.wave.seqera.io/library/bowtie2_htslib_samtools_pigz:e1c7af09e1513a6c' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('bowtie2')    , emit: index
    tuple val("${task.process}"), val('bowtie2'), eval("bowtie2 --version 2>&1 | sed -En 's/.*bowtie2-align-s(-v[0-9]+)? version //p'"), emit: versions_bowtie2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir bowtie2
    bowtie2-build $args --threads $task.cpus $fasta bowtie2/${fasta.baseName}
    """

    stub:
    """
    mkdir bowtie2
    touch bowtie2/${fasta.baseName}.{1..4}.bt2
    touch bowtie2/${fasta.baseName}.rev.{1,2}.bt2
    """
}
