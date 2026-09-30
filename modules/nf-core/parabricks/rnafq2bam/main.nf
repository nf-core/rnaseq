nextflow.enable.types = true

include { ReadsInput; StarLogs } from '../../types'

record ParabricksRnafq2bamResult {
    id:                 String
    meta:               Map
    raw_bams:           List<Path>
    bam_sorted:         Path?
    bam_sorted_aligned: Path?
    bam_unsorted:       Path?
    transcriptome_bam:  Path?
    unmapped:           List<Path>
    sam:                Path?
    junction:           Path?
    spl_junc_tab:       Path?
    read_per_gene_tab:  Path?
    wig:                List<Path>
    bedgraph:           List<Path>
    orig_bai:           Path?
    qc_metrics:         Path?
    duplicate_metrics:  Path?
    star:               StarLogs
}

process PARABRICKS_RNAFQ2BAM {
    tag "${sample.meta.id}"
    label 'process_high'
    label 'process_gpu'
    // needed by the module to work properly can be removed when fixed upstream - see: https://github.com/nf-core/modules/issues/7226
    stageInMode 'copy'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.7.1-1"

    input:
    sample: ReadsInput
    fasta: Path
    index: Path
    qc_metrics: Boolean
    mark_duplicates: Boolean

    output:
    record(
        id:                sample.id,
        meta:              sample.meta,
        raw_bams:          files("${prefix}.bam", optional: true).toSorted { f -> f.name },
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
        orig_bai:          file("${prefix}.bam.bai", optional: true),
        qc_metrics:        file("${prefix}_qc_metrics", optional: true),
        duplicate_metrics: file("${prefix}.duplicate-metrics.txt", optional: true),
        star:              record(
            log_final:    file("${prefix}.Log.final.out"),
            log_out:      file("${prefix}.Log.out"),
            log_progress: file("${prefix}.Log.progress.out"),
            tab:          files('*.tab', optional: true).toSorted { f -> f.name }
        )
    ) as ParabricksRnafq2bamResult

    topic:
    tuple(task.process, 'parabricks', eval("pbrun version 2>&1 | grep -Po '(?<=^pbrun: ).*'")) >> 'versions'

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${sample.meta.id}"

    def in_fq_command = sample.meta.single_end ? "--in-se-fq ${sample.reads.join(' ')}" : "--in-fq ${sample.reads.join(' ')}"
    def num_gpus = task.accelerator ? "--num-gpus ${task.accelerator.request}" : ''

    def qc_metrics_command = qc_metrics ? "--out-qc-metrics-dir ${prefix}_qc_metrics" : ""
    def duplicate_metrics_command = mark_duplicates ? "--out-duplicate-metrics ${prefix}.duplicate-metrics.txt" : "--no-markdups"

    """
    pbrun \\
        rna_fq2bam  \\
        --ref ${fasta} \\
        ${in_fq_command} \\
        --output-dir . \\
        --genome-lib-dir ${index} \\
        --out-bam ${prefix}.bam \\
        --logfile ${prefix}.Log.final.out \\
        --out-prefix ${prefix}. \\
        ${num_gpus} \\
        ${qc_metrics_command} \\
        ${duplicate_metrics_command} \\
        ${args}
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    prefix = task.ext.prefix ?: "${sample.meta.id}"
    def qc_metrics_output = qc_metrics ? "mkdir ${prefix}_qc_metrics" : ""
    def duplicate_metrics_output = mark_duplicates ? "touch ${prefix}.duplicate-metrics.txt" : ""
    """
    echo "" | gzip > ${prefix}.unmapped_1.fastq.gz
    echo "" | gzip > ${prefix}.unmapped_2.fastq.gz
    touch ${prefix}.bam
    touch ${prefix}.bam.bai
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
    ${qc_metrics_output}
    ${duplicate_metrics_output}
    """
}
