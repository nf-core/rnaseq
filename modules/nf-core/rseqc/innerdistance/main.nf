nextflow.enable.types = true

include { BamBaiInput; RseqcInnerDistance } from '../../types'

process RSEQC_INNERDISTANCE {
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
        id:     sample.id,
        meta:     sample.meta,
        distance: file("*distance.txt"),
        freq:     file("*freq.txt", optional: true),
        mean:     file("*mean.txt", optional: true),
        pdf:      file("*.pdf", optional: true),
        rscript:  file("*.r", optional: true)
    ) as RseqcInnerDistance

    topic:
    tuple(task.process, 'rseqc', eval('inner_distance.py --version | sed "s/inner_distance.py //"')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    if (!sample.meta.single_end) {
        """
        inner_distance.py \\
            -i $sample.bam \\
            -r $bed \\
            -o $prefix \\
            $args \\
            > stdout.txt
        head -n 2 stdout.txt > ${prefix}.inner_distance_mean.txt
        """
    } else {
        """
        echo "inner_distance.py doesn't support single-end data" > ${prefix}.inner_distance.txt
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    """
    touch ${prefix}.inner_distance.txt
    touch ${prefix}.inner_distance_freq.txt
    touch ${prefix}.inner_distance_mean.txt
    touch ${prefix}.inner_distance_plot.pdf
    touch ${prefix}.inner_distance_plot.r
    """
}
