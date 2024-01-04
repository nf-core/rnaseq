process RSEQC_JUNCTIONANNOTATION {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::rseqc=5.0.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rseqc:5.0.3--py39hf95cd2a_0' :
        'biocontainers/rseqc:5.0.3--py39hf95cd2a_0' }"

    input:
    tuple val(meta), path(bam)
    path  bed

    output:
    tuple val(meta), path("*.xls")         , emit: xls
    tuple val(meta), path("*.r")           , emit: rscript
    tuple val(meta), path("*.log")         , emit: log
    tuple val(meta), path("*.junction.bed"), optional:true, emit: bed
    tuple val(meta), path("*.Interact.bed"), optional:true, emit: interact_bed
    tuple val(meta), path("*junction.pdf") , optional:true, emit: pdf
    tuple val(meta), path("*events.pdf")   , optional:true, emit: events_pdf
    path  "versions.yml"                   , emit: versions

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
        2> ${prefix}.junction_annotation.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(junction_annotation.py --version | sed -e "s/junction_annotation.py //g")
    END_VERSIONS
    """
}
