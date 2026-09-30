nextflow.enable.types = true

include { BamBaiInput } from '../../types'

process RSEQC_READDISTRIBUTION {
    tag "$sample.meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6f/6f44b7933e2c2b1a340dc9485869974eb032d34e81af83716eb381964ee3e5e7/data' :
        'community.wave.seqera.io/library/rseqc_r-base:2e29d2dfda9cef15' }"

    input:
    sample: BamBaiInput
    bed: Path

    output:
    record(id: sample.id, meta: sample.meta, readdistribution: file("*.read_distribution.txt"))

    topic:
    tuple(task.process, 'rseqc', eval('read_distribution.py --version | sed "s/read_distribution.py //"')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    """
    read_distribution.py \\
        $args \\
        -i $sample.bam \\
        -r $bed \\
        > ${prefix}.read_distribution.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    """
    touch ${prefix}.read_distribution.txt
    """
}
