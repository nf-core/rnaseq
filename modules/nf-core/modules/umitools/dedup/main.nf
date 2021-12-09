process UMITOOLS_DEDUP {
    tag "$meta.id"
    label "process_medium"

    conda (params.enable_conda ? "bioconda::umi_tools=1.1.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/umi_tools:1.1.2--py38h4a8c8d9_0"
    } else {
        container "quay.io/biocontainers/umi_tools:1.1.2--py38h4a8c8d9_0"
    }

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path  "versions.yml"          , emit: versions

    script:
    def prefix   = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def args     = task.ext.args ?: ''
    def paired   = meta.single_end ? "" : "--paired"
    """
    umi_tools dedup \\
        -I $bam \\
        -S ${prefix}.bam \\
        $paired \\
        $args

    cat <<-END_VERSIONS > versions.yml
    ${task.process.tokenize(':').last()}:
        umitools: \$(umi_tools --version 2>&1 | sed 's/^.*UMI-tools version://; s/ *\$//')
    END_VERSIONS
    """
}
