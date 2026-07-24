process RSEQC_TIN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6f/6f44b7933e2c2b1a340dc9485869974eb032d34e81af83716eb381964ee3e5e7/data' :
        'community.wave.seqera.io/library/rseqc_r-base:2e29d2dfda9cef15' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  bed

    output:
    tuple val(meta), path("*.txt"), emit: txt
    tuple val(meta), path("*.xls"), emit: xls
    tuple val("${task.process}"), val('rseqc'), eval('tin.py --version | sed "s/tin.py //"'), emit: versions_rseqc, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    tin.py \\
        -i $bam \\
        -r $bed \\
        $args

    mv ${bam.baseName}.summary.txt ${prefix}.summary.txt
    mv ${bam.baseName}.tin.xls ${prefix}.tin.xls
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.summary.txt
    touch ${prefix}.tin.xls
    """
}
