process KALLISTO_QUANT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kallisto:0.48.0--h15996b6_2':
        'biocontainers/kallisto:0.48.0--h15996b6_2' }"

    input:
    tuple val(meta), path(reads)
    path index
    path gtf
    path chromosomes

    output:
    tuple val(meta), path("abundance.tsv"), emit: abundance
    tuple val(meta), path("abundance.h5") , emit: abundance_hdf5
    tuple val(meta), path("run_info.json"), emit: run_info
    tuple val(meta), path("*.log.txt")    , emit: log
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def single = meta.single_end ? "--single --fragment-length ${task.ext.fragment_len} --sd ${task.ext.sd}" : ""
    def gtf_input = gtf ? "--gtf ${gtf}" : ''
    def chromosomes_input = chromosomes ? "--chromosomes ${chromosomes}" : ''
    """
    kallisto quant \\
            --threads ${task.cpus} \\
            --index ${index} \\
            ${gtf_input} \\
            ${chromosomes_input} \\
            ${single} \\
            ${args} \\
            -o . \\
            ${reads} 2> >(tee -a ${prefix}.log.txt >&2)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallisto: \$(echo \$(kallisto version) | sed "s/kallisto, version //g" )
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallisto: \$(echo \$(kallisto version) | sed "s/kallisto, version //g" )
    END_VERSIONS
    """
}
