process RSEQC_JUNCTIONANNOTATION {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rseqc:5.0.3--py39hf95cd2a_0' :
        'biocontainers/rseqc:5.0.3--py39hf95cd2a_0' }"

    input:
    meta    : Map
    bam     : Path
    bed     : Path

    output:
    xls             : Path = file("*.xls")
    rscript         : Path = file("*.r")
    log             : Path = file("*.log")
    bed             : Path? = file("*.junction.bed")
    interact_bed    : Path? = file("*.Interact.bed")
    pdf             : Path? = file("*junction.pdf")
    events_pdf      : Path? = file("*events.pdf")

    topic:
    file('versions.yml') >> 'versions'
    file('*.log') >> 'logs'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    junction_annotation.py \\
        -i $bam \\
        -r $bed \\
        -o $prefix \\
        $args \\
        2> >(grep -v 'E::idx_find_and_load' | tee ${prefix}.junction_annotation.log >&2)
    """
}
