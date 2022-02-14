process RSEQC_JUNCTIONSATURATION {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::rseqc=3.0.1 'conda-forge::r-base>=3.5'" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rseqc:3.0.1--py37h516909a_1' :
        'quay.io/biocontainers/rseqc:3.0.1--py37h516909a_1' }"

    input:
    tuple val(meta), path(bam)
    path  bed

    output:
    tuple val(meta), path("*.pdf"), emit: pdf
    tuple val(meta), path("*.r")  , emit: rscript
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    junction_saturation.py \\
        -i $bam \\
        -r $bed \\
        -o $prefix \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(junction_saturation.py --version | sed -e "s/junction_saturation.py //g")
    END_VERSIONS
    """
}
