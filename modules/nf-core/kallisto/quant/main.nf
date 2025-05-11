process KALLISTO_QUANT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kallisto:0.48.0--h15996b6_2':
        'biocontainers/kallisto:0.48.0--h15996b6_2' }"

    input:
    meta                : Map
    reads               : List<Path>
    index               : Path
    gtf                 : Path
    chromosomes         : Path?
    fragment_length     : int
    fragment_length_sd  : int

    output:
    results     : Path = file("${prefix}")
    json_info   : Path = file("*.run_info.json")
    log         : Path = file("*.log")

    topic:
    file('versions.yml') >> 'versions'
    file('*.log') >> 'logs'

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
    """

    stub:
    """
    """
}
