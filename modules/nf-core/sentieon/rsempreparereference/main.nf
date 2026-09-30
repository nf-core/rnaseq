nextflow.enable.types = true

process SENTIEON_RSEMPREPAREREFERENCE {
    tag "$fasta"
    label 'process_high'
    label 'sentieon'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/39/39a3e1a85912520836ad054c8ac0497b463bb5170e0e907183dbd08509dad997/data' :
        'community.wave.seqera.io/library/rsem_sentieon:3e4315fa0b636313' }"

    input:
    fasta: Path
    gtf: Path

    stage:
    stageAs fasta, 'rsem/*'

    output:
    record(
        index:            file("rsem"),
        transcript_fasta: file("*transcripts.fa")
    )

    topic:
    tuple(task.process, 'rsem', eval('rsem-calculate-expression --version | sed -e "s/Current version: RSEM v//g"')) >> 'versions'
    tuple(task.process, 'star', eval('STAR --version | sed -e "s/STAR_//g"')) >> 'versions'
    tuple(task.process, 'sentieon', eval('sentieon driver --version 2>&1 | sed -e "s/sentieon-genomics-//g"')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args_list = args.tokenize()
    def args_no_star = args_list.findAll { arg -> !arg.contains('--star') }.join(' ')
    def star_symlink = """
        # Create symlink to sentieon in PATH for version detection
        ln -sf \$(which sentieon) ./STAR
        export PATH=".:\$PATH"
        """
    if (args_list.contains('--star')) {
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
            ${args_no_star} \\
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
