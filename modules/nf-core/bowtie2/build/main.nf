// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process BOWTIE2_BUILD {
    tag "$fasta"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d4/d4990a9e13c1aad32d60f19e0ced204f2f1e231fe856290ce75713505a1bfe44/data'
:         'community.wave.seqera.io/library/bowtie2_htslib_samtools_pigz:dbd18842ad7cb278' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('bowtie2')    , emit: index
    tuple val("${task.process}"), val('bowtie2'), eval("bowtie2 --version 2>&1 | sed -n 's/.*bowtie2-align-s version //p'"), emit: versions_bowtie2, topic: versions

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
