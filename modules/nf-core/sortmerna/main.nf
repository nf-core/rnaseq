process SORTMERNA {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/91/919d9c8f5f2c3221a94efe96b81bde0c953c13ebb0a1eca6b690b90666006cad/data' :
        'community.wave.seqera.io/library/sortmerna:4.3.7--b730cad73fc42b8e' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(fastas)
    tuple val(meta3), path(index)

    output:
    tuple val(meta), path("*non_rRNA.fastq.gz"), emit: reads, optional: true
    tuple val(meta), path("*.log")             , emit: log, optional: true
    tuple val(meta2), path("idx")              , emit: index, optional: true
    path  "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args  ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"

    def index_only    = args.contains('--index 1')? true : false
    def skip_index    = args.contains('--index 0')? true : false
    def paired_end    = reads instanceof List
    def paired_cmd    = ''
    def reads_args    = ''
    def out2_cmd      = ''
    def mv_cmd        = ''
    def reads_input   = ''
    def refs_input    = ''

    if (! index_only){
        reads_args = '--aligned rRNA_reads --fastx --other non_rRNA_reads'
        reads_input = paired_end ? reads.collect{"--reads $it"}.join(' ') : "--reads $reads"
        def n_fastq = paired_end ? reads.size() : 1
        if ( n_fastq == 1 ) {
            mv_cmd = """
            mv non_rRNA_reads.f*q.gz ${prefix}.non_rRNA.fastq.gz
            mv rRNA_reads.log ${prefix}.sortmerna.log
            """
        } else {
            mv_cmd = """
            mv non_rRNA_reads_fwd.f*q.gz ${prefix}_1.non_rRNA.fastq.gz
            mv non_rRNA_reads_rev.f*q.gz ${prefix}_2.non_rRNA.fastq.gz
            mv rRNA_reads.log ${prefix}.sortmerna.log
            """
            paired_cmd = "--paired_in"
            out2_cmd   = "--out2"
        }
    }
    """
    sortmerna \\
        ${'--ref '+fastas.join(' --ref ')} \\
        $refs_input \\
        $reads_input \\
        --threads $task.cpus \\
        --workdir . \\
        $reads_args \\
        $paired_cmd \\
        $out2_cmd \\
        $args

    $mv_cmd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sortmerna: \$(echo \$(sortmerna --version 2>&1) | sed 's/^.*SortMeRNA version //; s/ Build Date.*\$//')
    END_VERSIONS
    """

    stub:
    def args          = task.ext.args  ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"

    def index_only    = args.contains('--index 1')? true : false
    def paired_end    = reads instanceof List
    def paired_cmd    = ''
    def out2_cmd      = ''
    def mv_cmd        = ''
    def reads_input   = ''

    if (! index_only){
        reads_input = paired_end ? reads.collect{"--reads $it"}.join(' ') : "--reads $reads"
        def n_fastq = paired_end ? reads.size() : 1
        if ( n_fastq == 1 ) {
            mv_cmd = "touch ${prefix}.non_rRNA.fastq.gz"
        } else {
            mv_cmd = """
            touch ${prefix}_1.non_rRNA.fastq.gz
            touch ${prefix}_2.non_rRNA.fastq.gz
            """
        }
    }
    """
    $mv_cmd
    mkdir -p idx
    touch ${prefix}.sortmerna.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sortmerna: \$(echo \$(sortmerna --version 2>&1) | sed 's/^.*SortMeRNA version //; s/ Build Date.*\$//')
    END_VERSIONS
    """
}
