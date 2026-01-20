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
    tuple val(meta), path("${prefix}.Log.final.out")                    , emit: log_final
    path  "versions.yml"                                                , emit: versions

    tuple val(meta), path("${prefix}.bam")                              , optional:true, emit: bam
    tuple val(meta), path("${prefix}.bai")                              , optional:true, emit: bai
    tuple val(meta), path("${prefix}.sortedByCoord.out.bam")            , optional:true, emit: bam_sorted
    tuple val(meta), path("${prefix}.Aligned.sortedByCoord.out.bam")    , optional:true, emit: bam_sorted_aligned
    tuple val(meta), path('*toTranscriptome.out.bam')                   , optional:true, emit: bam_transcript
    tuple val(meta), path('*Aligned.unsort.out.bam')                    , optional:true, emit: bam_unsorted
    tuple val(meta), path('*fastq.gz')                                  , optional:true, emit: fastq
    tuple val(meta), path('*.tab')                                      , optional:true, emit: tab
    tuple val(meta), path('*.SJ.out.tab')                               , optional:true, emit: spl_junc_tab
    tuple val(meta), path('*.ReadsPerGene.out.tab')                     , optional:true, emit: read_per_gene_tab
    tuple val(meta), path('*.out.junction')                             , optional:true, emit: junction
    tuple val(meta), path('*.out.sam')                                  , optional:true, emit: sam
    tuple val(meta), path('*.wig')                                      , optional:true, emit: wig
    tuple val(meta), path('*.bg')                                       , optional:true, emit: bedgraph
    tuple val(meta), path("${prefix}_qc_metrics")                       , optional:true, emit: qc_metrics     
    tuple val(meta), path("${prefix}.duplicate-metrics.txt")            , optional:true, emit: duplicate_metrics

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
        ${num_gpus} \\
        ${qc_metrics_command} \\
        ${duplicate_metrics_command} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    prefix = task.ext.prefix ?: "${meta.id}"
    def qc_metrics_output = args.contains("--out-qc-metrics-dir") ? "mkdir ${prefix}_qc_metrics" : ""
    def duplicate_metrics_output = args.contains("--out-duplicate-metrics") ? "touch ${prefix}.duplicate-metrics.txt" : ""
    """
    echo "" | gzip > ${prefix}.unmapped_1.fastq.gz
    echo "" | gzip > ${prefix}.unmapped_2.fastq.gz
    touch ${prefix}.bam
    touch ${prefix}.bam.bai
    touch ${prefix}.Log.final.out
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

    # Capture the full version output once and store it in a variable
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """
}