nextflow.preview.types = true

record ReadDistributionResult {
    meta: Map
    txt:  Path
}

process RSEQC_READDISTRIBUTION {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6f/6f44b7933e2c2b1a340dc9485869974eb032d34e81af83716eb381964ee3e5e7/data' :
        'community.wave.seqera.io/library/rseqc_r-base:2e29d2dfda9cef15' }"

    input:
    (meta: Map, bam: Path, bai: Path): Record
    bed: Path

    output:
    record(meta: meta, txt: file("*.read_distribution.txt"))
    tuple val("${task.process}"), val('rseqc'), eval('read_distribution.py --version | sed "s/read_distribution.py //"'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    read_distribution.py \\
        $args \\
        -i $bam \\
        -r $bed \\
        > ${prefix}.read_distribution.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.read_distribution.txt
    """
}
