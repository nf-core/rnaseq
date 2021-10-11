process RSEQC_JUNCTIONSATURATION {
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
    tuple val(meta), path("*.pdf"), emit: pdf
    tuple val(meta), path("*.r")  , emit: rscript
    path  "versions.yml"          , emit: versions

    script:
    def prefix   = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def args     = task.ext.args ?: ''
    """
    junction_saturation.py \\
        -i $bam \\
        -r $bed \\
        -o $prefix \\
        $args

    cat <<-END_VERSIONS > versions.yml
    RSEQC_JUNCTIONSATURATION:
        rseqc: \$(junction_saturation.py --version | sed -e "s/junction_saturation.py //g")
    END_VERSIONS
    """
}
