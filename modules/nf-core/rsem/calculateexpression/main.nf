process RSEM_CALCULATEEXPRESSION {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-cf0123ef83b3c38c13e3b0696a3f285d3f20f15b:64aad4a4e144878400649e71f42105311be7ed87-0' :
        'biocontainers/mulled-v2-cf0123ef83b3c38c13e3b0696a3f285d3f20f15b:64aad4a4e144878400649e71f42105311be7ed87-0' }"

    input:
    meta    : Map
    reads   : List<Path>
    index   : Path

    output:
    counts_gene         : Path = file("*.genes.results")
    counts_transcript   : Path = file("*.isoforms.results")
    stat                : Path = file("*.stat")
    logs                : Path = file("*.log")

    bam_star        : Path? = file("*.STAR.genome.bam")
    bam_genome      : Path? = file("${prefix}.genome.bam")
    bam_transcript  : Path? = file("${prefix}.transcript.bam")

    topic:
    file('versions.yml') >> 'versions'
    file("*.stat") >> 'logs'
    file("*.log") >> 'logs'

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    def strandedness = ''
    if (meta.strandedness == 'forward') {
        strandedness = '--strandedness forward'
    } else if (meta.strandedness == 'reverse') {
        strandedness = '--strandedness reverse'
    }
    def paired_end = meta.single_end ? "" : "--paired-end"
    """
    INDEX=`find -L ./ -name "*.grp" | sed 's/\\.grp\$//'`
    rsem-calculate-expression \\
        --num-threads $task.cpus \\
        --temporary-folder ./tmp/ \\
        $strandedness \\
        $paired_end \\
        $args \\
        $reads \\
        \$INDEX \\
        $prefix
    """
}
