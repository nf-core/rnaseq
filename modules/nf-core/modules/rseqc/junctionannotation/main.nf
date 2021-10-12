process RSEQC_JUNCTIONANNOTATION {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::rseqc=3.0.1 'conda-forge::r-base>=3.5'" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/rseqc:3.0.1--py37h516909a_1"
    } else {
        container "quay.io/biocontainers/rseqc:3.0.1--py37h516909a_1"
    }

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

    script:
    def prefix   = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def args     = task.ext.args ?: ''
    """
    junction_annotation.py \\
        -i $bam \\
        -r $bed \\
        -o $prefix \\
        $args \\
        2> ${prefix}.junction_annotation.log

    cat <<-END_VERSIONS > versions.yml
    RSEQC_JUNCTIONANNOTATION:
        rseqc: \$(junction_annotation.py --version | sed -e "s/junction_annotation.py //g")
    END_VERSIONS
    """
}
