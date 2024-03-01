process SORTMERNA {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sortmerna:4.3.6--h9ee0642_0' :
        'biocontainers/sortmerna:4.3.6--h9ee0642_0' }"

    input:
    tuple val(meta), path(reads)
    path  fastas

    output:
    tuple val(meta), path("*non_rRNA.fastq.gz"), emit: reads
    tuple val(meta), path("*.log")     , emit: log
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args  ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def reads_input = reads instanceof List ? reads.collect{"--reads $it"}.join(' ') : "--reads $reads"
    def n_fastq     = reads instanceof List ? reads.size() : 1
    if ( n_fastq == 1 ) {
        mv_cmd = "mv non_rRNA_reads.f*q.gz ${prefix}.non_rRNA.fastq.gz"
        paired_cmd = ''
        out2_cmd   = ''
    } else {
        mv_cmd = """
        mv non_rRNA_reads_fwd.f*q.gz ${prefix}_1.non_rRNA.fastq.gz
        mv non_rRNA_reads_rev.f*q.gz ${prefix}_2.non_rRNA.fastq.gz
        """.stripIndent()
        paired_cmd = "--paired_in"
        out2_cmd   = "--out2"
    }
    """
    sortmerna \\
        ${'--ref '+fastas.join(' --ref ')} \\
        $reads_input \\
        --threads $task.cpus \\
        --workdir . \\
        --aligned rRNA_reads \\
        --fastx \\
        --other non_rRNA_reads \\
        $paired_cmd \\
        $out2_cmd \\
        $args

    $mv_cmd
    mv rRNA_reads.log ${prefix}.sortmerna.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sortmerna: \$(echo \$(sortmerna --version 2>&1) | sed 's/^.*SortMeRNA version //; s/ Build Date.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads_input = reads instanceof List ? reads.collect{"--reads $it"}.join(' ') : "--reads $reads"
    def n_fastq     = reads instanceof List ? reads.size() : 1
    if ( n_fastq == 1 ) {
        mv_cmd = "touch ${prefix}.non_rRNA.fastq.gz"
    } else {
        mv_cmd = """
        touch ${prefix}_1.non_rRNA.fastq.gz
        touch ${prefix}_2.non_rRNA.fastq.gz
        """.stripIndent()
    }
    """
    $mv_cmd
    touch ${prefix}.sortmerna.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sortmerna: \$(echo \$(sortmerna --version 2>&1) | sed 's/^.*SortMeRNA version //; s/ Build Date.*\$//')
    END_VERSIONS
    """
}
