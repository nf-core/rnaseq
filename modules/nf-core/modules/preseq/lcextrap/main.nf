process PRESEQ_LCEXTRAP {
    tag "$meta.id"
    label 'process_medium'
    label 'error_ignore'

    conda (params.enable_conda ? "bioconda::preseq=3.1.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/preseq:3.1.2--h06ef8b0_1"
    } else {
        container "quay.io/biocontainers/preseq:3.1.2--h06ef8b0_1"
    }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.ccurve.txt"), emit: ccurve
    tuple val(meta), path("*.log")       , emit: log
    path  "versions.yml"                 , emit: versions

    script:
    def prefix     = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def args       = task.ext.args ?: ''
    def paired_end = meta.single_end ? '' : '-pe'
    """
    preseq \\
        lc_extrap \\
        $args \\
        $paired_end \\
        -output ${prefix}.ccurve.txt \\
        $bam
    cp .command.err ${prefix}.command.log

    cat <<-END_VERSIONS > versions.yml
    ${task.process.tokenize(':').last()}:
        preseq: \$(echo \$(preseq 2>&1) | sed 's/^.*Version: //; s/Usage:.*\$//')
    END_VERSIONS
    """
}
