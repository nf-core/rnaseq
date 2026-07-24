process PARABRICKS_RNAFQ2BAM {
    tag "${meta.id}"
    label 'process_high'
    label 'process_gpu'
    // needed by the module to work properly can be removed when fixed upstream - see: https://github.com/nf-core/modules/issues/7226
    stageInMode 'copy'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1"

    input:
    tuple val(meta),  path(reads)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(index)
    val qc_metrics
    val mark_duplicates

    output:
    tuple val(meta), path("${prefix}.Log.final.out"),                           emit: log_final
    tuple val(meta), path("${prefix}.Log.out"),                                 emit: log_out
    tuple val(meta), path("${prefix}.Log.progress.out"),                        emit: log_progress
    tuple val(meta), path("${prefix}.bam"),                                     emit: bam,                  optional:true
    tuple val(meta), path("${prefix}.bam.bai"),                                 emit: bai,                  optional:true
    tuple val(meta), path("*.sortedByCoord.out.bam"),                           emit: bam_sorted,           optional:true
    tuple val(meta), path("*.Aligned.sortedByCoord.out.bam"),                   emit: bam_sorted_aligned,   optional:true
    tuple val(meta), path('*toTranscriptome.out.bam'),                          emit: bam_transcript,       optional:true
    tuple val(meta), path('*Aligned.unsort.out.bam'),                           emit: bam_unsorted,         optional:true
    tuple val(meta), path('*fastq.gz'),                                         emit: fastq,                optional:true
    tuple val(meta), path('*.tab'),                                             emit: tab,                  optional:true
    tuple val(meta), path('*.SJ.out.tab'),                                      emit: spl_junc_tab,         optional:true
    tuple val(meta), path('*.ReadsPerGene.out.tab'),                            emit: read_per_gene_tab,    optional:true
    tuple val(meta), path('*.out.junction'),                                    emit: junction,             optional:true
    tuple val(meta), path('*.out.sam'),                                         emit: sam,                  optional:true
    tuple val(meta), path('*.wig'),                                             emit: wig,                  optional:true
    tuple val(meta), path('*.bg'),                                              emit: bedgraph,             optional:true
    tuple val(meta), path("${prefix}_qc_metrics"),                              emit: qc_metrics,           optional:true
    tuple val(meta), path("${prefix}.duplicate-metrics.txt"),                   emit: duplicate_metrics,    optional:true
    tuple val("${task.process}"), val("parabricks"), eval("pbrun version 2>&1 | grep -Po '(?<=^pbrun: ).*'"),   emit: versions_parabricks,  topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def in_fq_command = meta.single_end ? "--in-se-fq ${reads}" : "--in-fq ${reads}"
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
    prefix = task.ext.prefix ?: "${meta.id}"
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
