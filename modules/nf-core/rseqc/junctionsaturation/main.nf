nextflow.enable.types = true

include { BamBaiInput; RseqcJunctionSaturation } from '../../types'

process RSEQC_JUNCTIONSATURATION {
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
    record(
        id:      sample.id,
        meta:    sample.meta,
        pdf:     file("*.pdf"),
        rscript: file("*.r")
    ) as RseqcJunctionSaturation

    topic:
    tuple(task.process, 'rseqc', eval('junction_saturation.py --version | sed "s/junction_saturation.py //"')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    """
    junction_saturation.py \\
        -i $sample.bam \\
        -r $bed \\
        -o $prefix \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    """
    touch ${prefix}.junctionSaturation_plot.pdf
    touch ${prefix}.junctionSaturation_plot.r
    """
}
