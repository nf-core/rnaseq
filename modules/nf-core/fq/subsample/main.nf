process FQ_SUBSAMPLE {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::fq=0.9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fq:0.9.1--h9ee0642_0':
        'biocontainers/fq:0.9.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    /* args requires:
        --probability <f64>: Probability read is kept, between 0 and 1. Mutually exclusive with record-count.
        --record-count <u64>: Number of records to keep. Mutually exclusive with probability
    */
    def args = task.ext.args ?: ''
    def prob_exists = args =~ /-p|--probability/
    def nrec_exists = args =~ /-n|--record-count/
    if ( !(prob_exists || nrec_exists) ){
        error "FQ/SUBSAMPLE requires --probability (-p) or --record-count (-n) specified in task.ext.args!"
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    def n_fastq = fastq instanceof List ? fastq.size() : 1
    log.debug "FQ/SUBSAMPLE found ${n_fastq} FASTQ files"
    if ( n_fastq == 1 ){
        fastq1_output = "--r1-dst ${prefix}.fastq.gz"
        fastq2_output = ""
    } else if ( n_fastq == 2 ){
        fastq1_output = "--r1-dst ${prefix}_R1.fastq.gz"
        fastq2_output = "--r2-dst ${prefix}_R2.fastq.gz"
    } else {
        error "FQ/SUBSAMPLE only accepts 1 or 2 FASTQ files!"
    }
    """
    fq subsample \\
        $args \\
        $fastq \\
        $fastq1_output \\
        $fastq2_output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fq: \$(echo \$(fq subsample --version | sed 's/fq-subsample //g'))
    END_VERSIONS
    """
}
