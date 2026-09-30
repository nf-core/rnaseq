nextflow.enable.types = true

process STAR_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/26/268b4c9c6cbf8fa6606c9b7fd4fafce18bf2c931d1a809a0ce51b105ec06c89d/data' :
        'community.wave.seqera.io/library/htslib_samtools_star_gawk:ae438e9a604351a4' }"

    input:
    record(id: String, meta: Map, reads: List<Path>)
    tuple(meta2: Map, index: Path)
    tuple(meta3: Map, gtf: Path?)
    star_ignore_sjdbgtf: Boolean

    stage:
    stageAs reads, 'input*/*'

    output:
    record(
        id:                id,
        meta:              meta,
        orig_bam:          files('*d.out.bam', optional: true).toSorted { f -> f.name },
        bam_sorted:        file("${prefix}.sortedByCoord.out.bam", optional: true),
        bam_sorted_aligned: file("${prefix}.Aligned.sortedByCoord.out.bam", optional: true),
        bam_unsorted:      file('*Aligned.unsort.out.bam', optional: true),
        transcriptome_bam: file('*toTranscriptome.out.bam', optional: true),
        unmapped:          files('*fastq.gz', optional: true).toSorted { f -> f.name },
        sam:               file('*.out.sam', optional: true),
        junction:          file('*.out.junction', optional: true),
        spl_junc_tab:      file('*.SJ.out.tab', optional: true),
        read_per_gene_tab: file('*.ReadsPerGene.out.tab', optional: true),
        wig:               files('*.wig', optional: true).toSorted { f -> f.name },
        bedgraph:          files('*.bg', optional: true).toSorted { f -> f.name },
        orig_bai:          null,
        qc_metrics:        null,
        duplicate_metrics: null,
        star:              record(
            log_final:    file('*Log.final.out'),
            log_out:      file('*Log.out'),
            log_progress: file('*Log.progress.out'),
            tab:          files('*.tab', optional: true).toSorted { f -> f.name }
        )
    )

    topic:
    tuple(task.process, 'star', eval('STAR --version | sed "s/STAR_//"')) >> 'versions'
    tuple(task.process, 'samtools', eval("samtools --version | sed -n '1s/samtools //p'")) >> 'versions'
    tuple(task.process, 'gawk', eval("gawk --version | sed -n '1s/GNU Awk \\([0-9.]*\\).*/\\1/p'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def read_pairs = reads.collate(2)
    def reads1 = meta.single_end ? reads : read_pairs.collect { pair -> pair[0] }.toList()
    def reads2 = meta.single_end ? [] : read_pairs.collect { pair -> pair[1] }.toList()
    def ignore_gtf      = star_ignore_sjdbgtf ? '' : "--sjdbGTFfile $gtf"
    attrRG          = args.contains("--outSAMattrRGline") ? "" : "--outSAMattrRGline 'ID:$prefix' 'SM:$prefix'"
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
