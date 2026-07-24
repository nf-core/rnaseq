process SENTIEON_RSEMPREPAREREFERENCE {
    tag "$fasta"
    label 'process_high'
    label 'sentieon'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/39/39a3e1a85912520836ad054c8ac0497b463bb5170e0e907183dbd08509dad997/data' :
        'community.wave.seqera.io/library/rsem_sentieon:3e4315fa0b636313' }"

    input:
    path fasta, stageAs: "rsem/*"
    path gtf

    output:
    path "rsem"           , emit: index
    path "*transcripts.fa", emit: transcript_fasta
    tuple val("${task.process}"), val('rsem'), eval('rsem-calculate-expression --version | sed -e "s/Current version: RSEM v//g"'), topic: versions, emit: versions_rsem
    tuple val("${task.process}"), val('star'), eval('STAR --version | sed -e "s/STAR_//g"'), topic: versions, emit: versions_star
    tuple val("${task.process}"), val('sentieon'), eval('sentieon driver --version 2>&1 | sed -e "s/sentieon-genomics-//g"'), topic: versions, emit: versions_sentieon

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args_list = args.tokenize()
    def star_symlink = """
        # Create symlink to sentieon in PATH for version detection
        ln -sf \$(which sentieon) ./STAR
        export PATH=".:\$PATH"
        """
    if (args_list.contains('--star')) {
        args_list.removeIf { arg -> arg.contains('--star') }
        def memory = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
        """
        $star_symlink

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
        """
    } else {
        """
        $star_symlink

        rsem-prepare-reference \\
            --gtf $gtf \\
            --num-threads $task.cpus \\
            $args \\
            $fasta \\
            rsem/genome

        cp rsem/genome.transcripts.fa .
        """
    }

    stub:
    """
    # Create symlink to sentieon in PATH for version detection
    ln -sf \$(which sentieon) ./STAR
    export PATH=".:\$PATH"

    mkdir -p rsem
    touch genome.transcripts.fa
    """
}
