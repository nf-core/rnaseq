// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process QUALIMAP_RNASEQ {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/43/43d53d976495e04b250e36cb8dd5dcbb0b1390542e9103071cd533a74269877d/data'
:         'community.wave.seqera.io/library/qualimap:2.3--e70af3f76915782f' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(gtf)

    output:
    tuple val(meta), path("${prefix}"), emit: results
    tuple val("${task.process}"), val('qualimap'), eval("qualimap -h | sed -n 's/^QualiMap v.//p'"), topic: versions, emit: versions_qualimap

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def paired_end = meta.single_end ? '' : '-pe'
    def memory = (task.memory.mega*0.8).intValue() + 'M'

    """
    unset DISPLAY
    mkdir -p tmp
    export _JAVA_OPTIONS=-Djava.io.tmpdir=./tmp
    qualimap \\
        --java-mem-size=${memory} \\
        rnaseq \\
        ${args} \\
        -bam ${bam} \\
        -gtf ${gtf} \\
        ${paired_end} \\
        -outdir ${prefix}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    """
}
