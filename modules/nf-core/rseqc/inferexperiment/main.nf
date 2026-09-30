nextflow.enable.types = true

include { BamBaiInput } from '../../types'

process RSEQC_INFEREXPERIMENT {
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
    record(id: sample.id, meta: sample.meta, inferexperiment: file("*.infer_experiment.txt"))

    topic:
    tuple(task.process, 'rseqc', eval('infer_experiment.py --version | sed "s/infer_experiment.py //"')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    """
    infer_experiment.py \\
        -i $sample.bam \\
        -r $bed \\
        $args \\
        > ${prefix}.infer_experiment.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    """
    touch ${prefix}.infer_experiment.txt
    """
}
