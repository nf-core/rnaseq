process SORTMERNA {
    tag "$meta.id"
    label "process_high"

    conda (params.enable_conda ? "bioconda::sortmerna=4.3.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sortmerna:4.3.4--h9ee0642_0' :
        'quay.io/biocontainers/sortmerna:4.3.4--h9ee0642_0' }"

    input:
    tuple val(meta), path(reads)
    path  fastas

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    tuple val(meta), path("*.log")     , emit: log
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        sortmerna \\
            ${'--ref '+fastas.join(' --ref ')} \\
            --reads $reads \\
            --threads $task.cpus \\
            --workdir . \\
            --aligned rRNA_reads \\
            --other non_rRNA_reads \\
            $args

        mv non_rRNA_reads.fq.gz ${prefix}.fastq.gz
        mv rRNA_reads.log ${prefix}.sortmerna.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sortmerna: \$(echo \$(sortmerna --version 2>&1) | sed 's/^.*SortMeRNA version //; s/ Build Date.*\$//')
        END_VERSIONS
        """
    } else {
        """
        sortmerna \\
            ${'--ref '+fastas.join(' --ref ')} \\
            --reads ${reads[0]} \\
            --reads ${reads[1]} \\
            --threads $task.cpus \\
            --workdir . \\
            --aligned rRNA_reads \\
            --other non_rRNA_reads \\
            --paired_in \\
            --out2 \\
            $args

        mv non_rRNA_reads_fwd.fq.gz ${prefix}_1.fastq.gz
        mv non_rRNA_reads_rev.fq.gz ${prefix}_2.fastq.gz
        mv rRNA_reads.log ${prefix}.sortmerna.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sortmerna: \$(echo \$(sortmerna --version 2>&1) | sed 's/^.*SortMeRNA version //; s/ Build Date.*\$//')
        END_VERSIONS
        """
    }
}
