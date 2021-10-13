process RSEM_PREPAREREFERENCE {
    tag "$fasta"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::rsem=1.3.3 bioconda::star=2.7.6a" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-cf0123ef83b3c38c13e3b0696a3f285d3f20f15b:606b713ec440e799d53a2b51a6e79dbfd28ecf3e-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-cf0123ef83b3c38c13e3b0696a3f285d3f20f15b:606b713ec440e799d53a2b51a6e79dbfd28ecf3e-0"
    }

    input:
    path fasta, stageAs: "rsem/*"
    path gtf

    output:
    path "rsem"                , emit: index
    path "rsem/*transcripts.fa", emit: transcript_fasta
    path "versions.yml"        , emit: versions

    script:
    def args     = (task.ext.args ?: '').tokenize()
    def args2    = task.ext.args2 ?: ''
    if (args.contains('--star')) {
        args.removeIf { it.contains('--star') }
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
            ${args.join(' ')} \\
            $fasta \\
            rsem/genome

        cat <<-END_VERSIONS > versions.yml
        ${task.process.tokenize(':').last()}:
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

        cat <<-END_VERSIONS > versions.yml
        ${task.process.tokenize(':').last()}:
            rsem: \$(rsem-calculate-expression --version | sed -e "s/Current version: RSEM v//g")
            star: \$(STAR --version | sed -e "s/STAR_//g")
        END_VERSIONS
        """
    }
}
