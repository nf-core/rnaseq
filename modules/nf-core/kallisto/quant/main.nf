nextflow.enable.types = true


record KallistoQuantInput {
    id:    String
    meta:  Map
    reads: List<Path>
    args:  String?
}

record KallistoQuantSample {
    id:        String
    meta:      Map
    quant_dir: Path
    json_info: Path
    log:       Path
}

process KALLISTO_QUANT {
    tag "$sample.meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e2/e21d2cff2526b0995996977c057f0c17844073781f86a18b15fef178a97ed7cb/data':
        'community.wave.seqera.io/library/kallisto:0.52.0--31c771060d82d25c' }"

    input:
    sample: KallistoQuantInput
    index: Path
    gtf: Path?
    chromosomes: Path?
    fragment_length: Integer?
    fragment_length_sd: Integer?

    output:
    record(
        id:        sample.id,
        meta:      sample.meta,
        quant_dir: file("${prefix}"),
        json_info: file("*.run_info.json"),
        log:       file("*.log")
    ) as KallistoQuantSample

    topic:
    tuple(task.process, 'kallisto', eval("kallisto version | sed 's/.*version //'")) >> 'versions'

    script:
    def args = task.ext.args ?: sample.args ?: ''
    prefix = task.ext.prefix ?: "${sample.meta.id}"
    def gtf_input = gtf ? "--gtf ${gtf}" : ''
    def chromosomes_input = chromosomes ? "--chromosomes ${chromosomes}" : ''

    def single_end_params = ''
    if (sample.meta.single_end) {
        if (!("${fragment_length}" =~ /^\d+$/)) {
            error "fragment_length must be set and numeric for single-end data"
        }
        if (!("${fragment_length_sd}" =~ /^\d+$/)) {
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
            ${sample.reads.join(' ')} 2>| >(tee -a ${prefix}/kallisto_quant.log >&2)

    cp ${prefix}/kallisto_quant.log ${prefix}.log
    cp ${prefix}/run_info.json ${prefix}.run_info.json
    """

    stub:
    prefix = task.ext.prefix ?: "${sample.meta.id}"

    """
    mkdir -p $prefix
    touch ${prefix}.log
    touch ${prefix}.run_info.json
    """
}
