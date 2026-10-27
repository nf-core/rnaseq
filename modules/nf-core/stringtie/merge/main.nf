process STRINGTIE_MERGE {
    tag "${meta.id}"
    label 'process_medium'

    // Note: 2.7X indices incompatible with AWS iGenomes.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3f/3facd74a0f728c9bb9e9a731b58c343895d2dbdfeb812ce5747f701103fc61cf/data' :
        'community.wave.seqera.io/library/stringtie:3.0.3--e8043d00caecd051' }"

    input:
    tuple val(meta), path(gtf)
    tuple val(meta2), path(annotation_gtf)

    output:
    tuple val(meta), path("${prefix}.gtf"), emit: merged_gtf
    tuple val("${task.process}"), val('stringtie'), eval('stringtie --version'), emit: versions_stringtie, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def reference = annotation_gtf ? "-G ${annotation_gtf}" : ""
    """
    stringtie \\
        --merge \\
        ${gtf} \\
        ${reference} \\
        -o ${prefix}.gtf \\
        -p ${task.cpus} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gtf
    """
}
