process SAMTOOLS_FASTQ {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta), path(input)
    val(interleave)

    output:
    tuple val(meta), path("*_{1,2}.fastq.gz")      , optional:true, emit: fastq
    tuple val(meta), path("*_interleaved.fastq")   , optional:true, emit: interleaved
    tuple val(meta), path("*_singleton.fastq.gz")  , optional:true, emit: singleton
    tuple val(meta), path("*_other.fastq.gz")      , optional:true, emit: other
    path  "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output = ( interleave && ! meta.single_end ) ? "> ${prefix}_interleaved.fastq" :
        meta.single_end ? "-1 ${prefix}_1.fastq.gz -s ${prefix}_singleton.fastq.gz" :
        "-1 ${prefix}_1.fastq.gz -2 ${prefix}_2.fastq.gz -s ${prefix}_singleton.fastq.gz"
    """
    # Note: --threads value represents *additional* CPUs to allocate (total CPUs = 1 + --threads).
    samtools \\
        fastq \\
        $args \\
        --threads ${task.cpus-1} \\
        -0 ${prefix}_other.fastq.gz \\
        $input \\
        $output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output = ( interleave && ! meta.single_end ) ? "touch ${prefix}_interleaved.fastq" :
        meta.single_end ? "echo | gzip > ${prefix}_1.fastq.gz && echo | gzip > ${prefix}_singleton.fastq.gz" :
        "echo | gzip > ${prefix}_1.fastq.gz && echo | gzip > ${prefix}_2.fastq.gz && echo | gzip > ${prefix}_singleton.fastq.gz"
    """
    ${output}
    echo | gzip > ${prefix}_other.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
