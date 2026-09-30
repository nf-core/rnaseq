nextflow.enable.types = true

include { FastaGtfInput } from '../../types'

process RSEM_PREPAREREFERENCE {
    tag "$sample.fasta"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/23/23651ffd6a171ef3ba867cb97ef615f6dd6be39158df9466fe92b5e844cd7d59/data' :
        'community.wave.seqera.io/library/rsem_star:5acb4e8c03239c32' }"

    input:
    sample: FastaGtfInput

    stage:
    stageAs sample.fasta, 'rsem/*'

    output:
    record(
        id:               sample.id,
        meta:             sample.meta,
        index:            file("rsem"),
        transcript_fasta: file("*transcripts.fa")
    )

    topic:
    tuple(task.process, 'rsem', eval('rsem-calculate-expression --version | sed -e "s/Current version: RSEM v//g"')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args_list = args.tokenize()
    def args_no_star = args_list.findAll { arg -> !arg.contains('--star') }.join(' ')
    if (args_list.contains('--star')) {
        def memory = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
        """
        STAR \\
            --runMode genomeGenerate \\
            --genomeDir rsem/ \\
            --genomeFastaFiles $sample.fasta \\
            --sjdbGTFfile $sample.gtf \\
            --runThreadN $task.cpus \\
            $memory \\
            $args2

        rsem-prepare-reference \\
            --gtf $sample.gtf \\
            --num-threads $task.cpus \\
            ${args_no_star} \\
            $sample.fasta \\
            rsem/genome

        cp rsem/genome.transcripts.fa .
        """
    } else {
        """
        rsem-prepare-reference \\
            --gtf $sample.gtf \\
            --num-threads $task.cpus \\
            $args \\
            $sample.fasta \\
            rsem/genome

        cp rsem/genome.transcripts.fa .
        """
    }

    stub:
    """
    touch genome.transcripts.fa
    """
}
