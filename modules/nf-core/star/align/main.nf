nextflow.preview.types = true

record StarAlignResult {
    meta:               Map
    bam:                Path?
    bam_sorted:         Path?
    bam_sorted_aligned: Path?
    bam_transcript:     Path?
    bam_unsorted:       Path?
    log_final:          Path
    log_out:            Path
    log_progress:       Path
    fastq:              Path?
    tab:                Path?
    spl_junc_tab:       Path?
    read_per_gene_tab:  Path?
    junction:           Path?
    sam:                Path?
    wig:                Path?
    bedgraph:           Path?
}

process STAR_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/26/268b4c9c6cbf8fa6606c9b7fd4fafce18bf2c931d1a809a0ce51b105ec06c89d/data' :
        'community.wave.seqera.io/library/htslib_samtools_star_gawk:ae438e9a604351a4' }"

    input:
    (meta: Map, reads: Path): Record
    (meta2: Map, index: Path): Record

    stage:
    stageAs(reads, 'input*/*')
    (meta3: Map, gtf: Path): Record
    star_ignore_sjdbgtf: String?
    seq_platform: String?
    seq_center: String?

    output:
    record(
        meta:               meta,
        bam:                file('*d.out.bam', optional: true),
        bam_sorted:         file("${prefix}.sortedByCoord.out.bam", optional: true),
        bam_sorted_aligned: file("${prefix}.Aligned.sortedByCoord.out.bam", optional: true),
        bam_transcript:     file('*toTranscriptome.out.bam', optional: true),
        bam_unsorted:       file('*Aligned.unsort.out.bam', optional: true),
        log_final:          file('*Log.final.out'),
        log_out:            file('*Log.out'),
        log_progress:       file('*Log.progress.out'),
        fastq:              file('*fastq.gz', optional: true),
        tab:                file('*.tab', optional: true),
        spl_junc_tab:       file('*.SJ.out.tab', optional: true),
        read_per_gene_tab:  file('*.ReadsPerGene.out.tab', optional: true),
        junction:           file('*.out.junction', optional: true),
        sam:                file('*.out.sam', optional: true),
        wig:                file('*.wig', optional: true),
        bedgraph:           file('*.bg', optional: true)
    )
    tuple val("${task.process}"), val('star'), eval('STAR --version | sed -e "s/STAR_//g"'), topic: versions
    tuple val("${task.process}"), val('samtools'), eval('echo $(samtools --version 2>&1) | sed "s/^.*samtools //; s/Using.*$//"'), topic: versions
    tuple val("${task.process}"), val('gawk'), eval('echo $(gawk --version 2>&1) | sed "s/^.*GNU Awk //; s/, .*$//"'), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def reads1 = []
    def reads2 = []
    meta.single_end ? [reads].flatten().each{ read -> reads1 << read} : reads.eachWithIndex{ v, ix -> ( ix & 1 ? reads2 : reads1) << v }
    def ignore_gtf      = star_ignore_sjdbgtf ? '' : "--sjdbGTFfile $gtf"
    def seq_platform_arg  = seq_platform ? "'PL:$seq_platform'" : ""
    def seq_center_arg    = seq_center ? "'CN:$seq_center'" : ""
    attrRG          = args.contains("--outSAMattrRGline") ? "" : "--outSAMattrRGline 'ID:$prefix' $seq_center_arg 'SM:$prefix' $seq_platform_arg"
    def out_sam_type    = (args.contains('--outSAMtype')) ? '' : '--outSAMtype BAM Unsorted'
    mv_unsorted_bam = (args.contains('--outSAMtype BAM Unsorted SortedByCoordinate')) ? "mv ${prefix}.Aligned.out.bam ${prefix}.Aligned.unsort.out.bam" : ''
    """
    STAR \\
        --genomeDir $index \\
        --readFilesIn ${reads1.join(",")} ${reads2.join(",")} \\
        --runThreadN $task.cpus \\
        --outFileNamePrefix $prefix. \\
        $out_sam_type \\
        $ignore_gtf \\
        $attrRG \\
        $args

    $mv_unsorted_bam

    if [ -f ${prefix}.Unmapped.out.mate1 ]; then
        mv ${prefix}.Unmapped.out.mate1 ${prefix}.unmapped_1.fastq
        gzip ${prefix}.unmapped_1.fastq
    fi
    if [ -f ${prefix}.Unmapped.out.mate2 ]; then
        mv ${prefix}.Unmapped.out.mate2 ${prefix}.unmapped_2.fastq
        gzip ${prefix}.unmapped_2.fastq
    fi

    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.unmapped_1.fastq.gz
    echo "" | gzip > ${prefix}.unmapped_2.fastq.gz
    touch ${prefix}Xd.out.bam
    touch ${prefix}.Log.final.out
    touch ${prefix}.Log.out
    touch ${prefix}.Log.progress.out
    touch ${prefix}.sortedByCoord.out.bam
    touch ${prefix}.toTranscriptome.out.bam
    touch ${prefix}.Aligned.unsort.out.bam
    touch ${prefix}.Aligned.sortedByCoord.out.bam
    touch ${prefix}.tab
    touch ${prefix}.SJ.out.tab
    touch ${prefix}.ReadsPerGene.out.tab
    touch ${prefix}.Chimeric.out.junction
    touch ${prefix}.out.sam
    touch ${prefix}.Signal.UniqueMultiple.str1.out.wig
    touch ${prefix}.Signal.UniqueMultiple.str1.out.bg
    """
}
