nextflow.preview.types = true

record ReadDuplicationResult {
    meta:    Map
    seq_xls: Path
    pos_xls: Path
    pdf:     Path
    rscript: Path
}

process RSEQC_READDUPLICATION {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6f/6f44b7933e2c2b1a340dc9485869974eb032d34e81af83716eb381964ee3e5e7/data' :
        'community.wave.seqera.io/library/rseqc_r-base:2e29d2dfda9cef15' }"

    input:
    (meta: Map, bam: Path, bai: Path): Record

    output:
    record(
        meta:    meta,
        seq_xls: file("*seq.DupRate.xls"),
        pos_xls: file("*pos.DupRate.xls"),
        pdf:     file("*.pdf"),
        rscript: file("*.r")
    )
    tuple val("${task.process}"), val('rseqc'), eval('read_duplication.py --version | sed "s/read_duplication.py //"'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    read_duplication.py \\
        -i $bam \\
        -o $prefix \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.seq.DupRate.xls
    touch ${prefix}.pos.DupRate.xls
    touch ${prefix}.DupRate_plot.pdf
    touch ${prefix}.DupRate_plot.r
    """
}
