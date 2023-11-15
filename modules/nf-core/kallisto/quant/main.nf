process KALLISTO_QUANT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kallisto:0.48.0--h15996b6_2':
        'biocontainers/kallisto:0.48.0--h15996b6_2' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    path gtf
    path chromosomes
    val fragment_length
    val fragment_length_sd

    output:
    tuple val(meta), path("${prefix}")        , emit: results
    tuple val(meta), path("*.run_info.json")  , emit: json_info
    tuple val(meta), path("*.log")            , emit: log
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def gtf_input = gtf ? "--gtf ${gtf}" : ''
    def chromosomes_input = chromosomes ? "--chromosomes ${chromosomes}" : ''

    def single_end_params = ''
    if (meta.single_end) {
        if (!(fragment_length =~ /^\d+$/)) {
            error "fragment_length must be set and numeric for single-end data"
        }
        if (!(fragment_length_sd =~ /^\d+$/)) {
            error "fragment_length_sd must be set and numeric for single-end data"
        }
        single_end_params = "--single --fragment-length=${fragment_length} --sd=${fragment_length_sd}"
    }

    def strandedness = ''
    if (!args.contains('--fr-stranded') && !args.contains('--rf-stranded')) {
        strandedness =  (meta.strandedness == 'forward') ? '--fr-stranded' :
                        (meta.strandedness == 'reverse') ? '--rf-stranded' : ''
    }

    """
    mkdir -p $prefix && kallisto quant \\
            --threads ${task.cpus} \\
            --index ${index} \\
            ${gtf_input} \\
            ${chromosomes_input} \\
            ${single_end_params} \\
            ${strandedness} \\
            ${args} \\
            -o $prefix \\
            ${reads} 2> >(tee -a ${prefix}/kallisto_quant.log >&2)

    cp ${prefix}/kallisto_quant.log ${prefix}.log
    cp ${prefix}/run_info.json ${prefix}.run_info.json

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
