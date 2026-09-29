nextflow.enable.types = true

process RSEQC_JUNCTIONANNOTATION {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6f/6f44b7933e2c2b1a340dc9485869974eb032d34e81af83716eb381964ee3e5e7/data' :
        'community.wave.seqera.io/library/rseqc_r-base:2e29d2dfda9cef15' }"

    input:
    tuple(meta: Map, bam: Path, bai: Path)
    bed: Path

    output:
    record(
        meta:         meta,
        bed:          file("*.junction.bed", optional: true),
        interact_bed: file("*.Interact.bed", optional: true),
        xls:          file("*.xls"),
        pdf:          file("*junction.pdf", optional: true),
        events_pdf:   file("*events.pdf", optional: true),
        rscript:      file("*.r"),
        log:          file("*.log")
    )

    topic:
    tuple(task.process, 'rseqc', eval('junction_annotation.py --version | sed "s/junction_annotation.py //"')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    junction_annotation.py \\
        -i $bam \\
        -r $bed \\
        -o $prefix \\
        $args \\
        2>| >(grep -v 'E::idx_find_and_load' | tee ${prefix}.junction_annotation.log >&2)
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.junction.xls
    touch ${prefix}.junction_plot.r
    touch ${prefix}.junction_annotation.log
    touch ${prefix}.junction.bed
    touch ${prefix}.Interact.bed
    touch ${prefix}.junction.pdf
    touch ${prefix}.events.pdf
    """
}
