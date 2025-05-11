process SORTMERNA {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sortmerna:4.3.6--h9ee0642_0' :
        'biocontainers/sortmerna:4.3.6--h9ee0642_0' }"

    input:
    meta    : Map
    reads   : List<Path>
    fastas  : List<Path>
    index   : Path?

    output:
    reads   : List<Path> = files("*non_rRNA.fastq.gz")
    log     : Path? = file("*.log")
    index   : Path? = file("idx")

    topic:
    file('versions.yml') >> 'versions'
    file('*.log') >> 'logs'

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
    """
}
