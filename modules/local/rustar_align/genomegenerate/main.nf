process RUSTAR_GENOMEGENERATE {
    tag "$fasta"
    label 'process_high'

    container "ghcr.io/scverse/rustar-aligner:dev"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(gtf)

    output:
    tuple val(meta), path("star")  , emit: index
    tuple val("${task.process}"), val('rustar-aligner'), eval("rustar-aligner --version | sed -n '1{s/^rustar-aligner //;p}'"), emit: versions_rustar, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def args_list   = args.tokenize()
    def include_gtf = gtf ? "--sjdbGTFfile $gtf" : ''
    def memory      = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    // Heuristic mirrors STAR_GENOMEGENERATE's gawk/samtools-faidx pipeline, but
    // computed in Groovy so we don't need samtools+gawk in the rustar container.
    // Approximating genome length with the on-disk fasta size is within 1-2% of
    // the true base count and is well inside the floor() rounding of log2(len)/2-1.
    def auto_sa_index = ''
    if (!args_list.contains('--genomeSAindexNbases')) {
        def genome_size = fasta.size()
        def computed = Math.floor(Math.log(genome_size as double) / Math.log(2) / 2 - 1) as int
        def num_bases = Math.min(14, Math.max(1, computed))
        auto_sa_index = "--genomeSAindexNbases ${num_bases}"
    }
    """
    mkdir star
    rustar-aligner \\
        --runMode genomeGenerate \\
        --genomeDir star/ \\
        --genomeFastaFiles $fasta \\
        $include_gtf \\
        --runThreadN $task.cpus \\
        $auto_sa_index \\
        $memory \\
        $args
    """

    stub:
    if (gtf) {
        """
        mkdir star
        touch star/Genome
        touch star/Log.out
        touch star/SA
        touch star/SAindex
        touch star/chrLength.txt
        touch star/chrName.txt
        touch star/chrNameLength.txt
        touch star/chrStart.txt
        touch star/exonGeTrInfo.tab
        touch star/exonInfo.tab
        touch star/geneInfo.tab
        touch star/genomeParameters.txt
        touch star/sjdbInfo.txt
        touch star/sjdbList.fromGTF.out.tab
        touch star/sjdbList.out.tab
        touch star/transcriptInfo.tab
        """
    } else {
        """
        mkdir star
        touch star/Genome
        touch star/Log.out
        touch star/SA
        touch star/SAindex
        touch star/chrLength.txt
        touch star/chrName.txt
        touch star/chrNameLength.txt
        touch star/chrStart.txt
        touch star/genomeParameters.txt
        """
    }
}
