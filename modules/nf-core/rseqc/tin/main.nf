nextflow.enable.types = true

include { BamBaiInput; RseqcTin } from '../../types'

process RSEQC_TIN {
    tag "$sample.meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6f/6f44b7933e2c2b1a340dc9485869974eb032d34e81af83716eb381964ee3e5e7/data' :
        'community.wave.seqera.io/library/rseqc_r-base:2e29d2dfda9cef15' }"

    input:
    sample: BamBaiInput
    bed: Path

    output:
    record(
        id:   sample.id,
        meta: sample.meta,
        txt:  file("*.txt"),
        xls:  file("*.xls")
    ) as RseqcTin

    topic:
    tuple(task.process, 'rseqc', eval('tin.py --version | sed "s/tin.py //"')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    """
    tin.py \\
        -i $sample.bam \\
        -r $bed \\
        $args

    mv ${sample.bam.baseName}.summary.txt ${prefix}.summary.txt
    mv ${sample.bam.baseName}.tin.xls ${prefix}.tin.xls
    """

    stub:
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    """
    touch ${prefix}.summary.txt
    touch ${prefix}.tin.xls
    """
}
