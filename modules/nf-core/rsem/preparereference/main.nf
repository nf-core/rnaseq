process RSEM_PREPAREREFERENCE {
    tag "$fasta"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-cf0123ef83b3c38c13e3b0696a3f285d3f20f15b:64aad4a4e144878400649e71f42105311be7ed87-0' :
        'biocontainers/mulled-v2-cf0123ef83b3c38c13e3b0696a3f285d3f20f15b:64aad4a4e144878400649e71f42105311be7ed87-0' }"

    input:
    path fasta, stageAs: "rsem/*"
    path gtf

    output:
    path "rsem"           , emit: index
    path "*transcripts.fa", emit: transcript_fasta
    path "versions.yml"   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args_list = args.tokenize()
    if (args_list.contains('--star')) {
        args_list.removeIf { it.contains('--star') }
        def memory = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
        """
        STAR \\
            --runMode genomeGenerate \\
            --genomeDir rsem/ \\
            --genomeFastaFiles $fasta \\
            --sjdbGTFfile $gtf \\
            --runThreadN $task.cpus \\
            $memory \\
            $args2

        rsem-prepare-reference \\
            --gtf $gtf \\
            --num-threads $task.cpus \\
            ${args_list.join(' ')} \\
            $fasta \\
            rsem/genome

        cp rsem/genome.transcripts.fa .

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            rsem: \$(rsem-calculate-expression --version | sed -e "s/Current version: RSEM v//g")
            star: \$(STAR --version | sed -e "s/STAR_//g")
        END_VERSIONS
        """
    } else {
        """
        rsem-prepare-reference \\
            --gtf $gtf \\
            --num-threads $task.cpus \\
            $args \\
            $fasta \\
            rsem/genome

        cp rsem/genome.transcripts.fa .

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            rsem: \$(rsem-calculate-expression --version | sed -e "s/Current version: RSEM v//g")
            star: \$(STAR --version | sed -e "s/STAR_//g")
        END_VERSIONS
        """
    }
}
