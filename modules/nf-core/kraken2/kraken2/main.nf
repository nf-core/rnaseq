nextflow.enable.types = true

include { ReadsInput; Kraken2Result } from '../../types'

process KRAKEN2_KRAKEN2 {
    tag "$sample.meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f827dcea51be6b5c32255167caa2dfb65607caecdc8b067abd6b71c267e2e82/data' :
        'community.wave.seqera.io/library/kraken2_coreutils_pigz:920ecc6b96e2ba71' }"

    input:
    sample: ReadsInput
    db: Path
    save_output_fastqs: Boolean
    save_reads_assignment: Boolean

    output:
    record(
        id:                          sample.id,
        meta:                        sample.meta,
        report:                      file('*report.txt'),
        classified_reads_fastq:      files('*.classified{.,_}*', optional: true),
        unclassified_reads_fastq:    files('*.unclassified{.,_}*', optional: true),
        classified_reads_assignment: file('*classifiedreads.txt', optional: true)
    ) as Kraken2Result

    topic:
    tuple(task.process, 'kraken2', eval('kraken2 --version 2>&1 | head -1 | sed "s/^.*Kraken version //; s/ .*//"')) >> 'versions'
    tuple(task.process, 'pigz', eval('pigz --version 2>&1 | sed "s/pigz //g"')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    def paired       = sample.meta.single_end ? "" : "--paired"
    def classified   = sample.meta.single_end ? "${prefix}.classified.fastq"   : "${prefix}.classified#.fastq"
    def unclassified = sample.meta.single_end ? "${prefix}.unclassified.fastq" : "${prefix}.unclassified#.fastq"
    def classified_option = save_output_fastqs ? "--classified-out ${classified}" : ""
    def unclassified_option = save_output_fastqs ? "--unclassified-out ${unclassified}" : ""
    def readclassification_option = save_reads_assignment ? "--output ${prefix}.kraken2.classifiedreads.txt" : "--output /dev/null"
    def compress_reads_command = save_output_fastqs ? "pigz -p $task.cpus *.fastq" : ""

    """
    kraken2 \\
        --db $db \\
        --threads $task.cpus \\
        --report ${prefix}.kraken2.report.txt \\
        --gzip-compressed \\
        $unclassified_option \\
        $classified_option \\
        $readclassification_option \\
        $paired \\
        $args \\
        ${sample.reads.join(' ')}

    $compress_reads_command
    """

    stub:
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    def classified   = sample.meta.single_end ? "${prefix}.classified.fastq.gz"   : "${prefix}.classified_1.fastq.gz ${prefix}.classified_2.fastq.gz"
    def unclassified = sample.meta.single_end ? "${prefix}.unclassified.fastq.gz" : "${prefix}.unclassified_1.fastq.gz ${prefix}.unclassified_2.fastq.gz"

    """
    touch ${prefix}.kraken2.report.txt
    if [ "$save_output_fastqs" == "true" ]; then
        touch $classified
        touch $unclassified
    fi
    if [ "$save_reads_assignment" == "true" ]; then
        touch ${prefix}.kraken2.classifiedreads.txt
    fi
    """

}
