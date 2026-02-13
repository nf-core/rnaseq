nextflow.preview.types = true

record JunctionAnnotationResult {
    meta:         Map
    bed:          Path?
    interact_bed: Path?
    xls:          Path
    log:          Path
    pdf:          Path?
    events_pdf:   Path?
    rscript:      Path
}

process RSEQC_JUNCTIONANNOTATION {
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
    record(
        meta:         meta,
        bed:          file("*.junction.bed", optional: true),
        interact_bed: file("*.Interact.bed", optional: true),
        xls:          file("*.xls"),
        log:          file("*.log"),
        pdf:          file("*junction.pdf",  optional: true),
        events_pdf:   file("*events.pdf",    optional: true),
        rscript:      file("*.r")
    )
    tuple val("${task.process}"), val('rseqc'), eval('junction_annotation.py --version | sed "s/junction_annotation.py //"'), topic: versions

    when:
    task.ext.when == null || task.ext.when

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
