process KALLISTO_QUANT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e2/e21d2cff2526b0995996977c057f0c17844073781f86a18b15fef178a97ed7cb/data':
        'community.wave.seqera.io/library/kallisto:0.52.0--31c771060d82d25c' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index), path(gtf), path(chromosomes)
    val fragment_length
    val fragment_length_sd

    output:
    tuple val(meta), path("${prefix}")        , emit: results
    tuple val(meta), path("*.run_info.json")  , emit: json_info
    tuple val(meta), path("*.log")            , emit: log
    tuple val("${task.process}"), val('kallisto'), eval("kallisto version | sed 's/.*version //'"), emit: versions_kallisto, topic: versions

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

    """
    mkdir -p $prefix && kallisto quant \\
            --threads ${task.cpus} \\
            --index ${index} \\
            ${gtf_input} \\
            ${chromosomes_input} \\
            ${single_end_params} \\
            ${args} \\
            -o $prefix \\
            ${reads} 2>| >(tee -a ${prefix}/kallisto_quant.log >&2)

    cp ${prefix}/kallisto_quant.log ${prefix}.log
    cp ${prefix}/run_info.json ${prefix}.run_info.json
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p $prefix
    touch ${prefix}.log
    touch ${prefix}.run_info.json
    """
}
