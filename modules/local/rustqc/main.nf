process RUSTQC {
    tag "$meta.id"
    label 'process_medium'

    container "ghcr.io/ewels/rustqc:dev"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(gtf)

    output:
    // dupRadar outputs
    tuple val(meta), path("*_duprateExpDens.png")          , emit: scatter2d         , optional: true
    tuple val(meta), path("*_duprateExpBoxplot.png")       , emit: boxplot           , optional: true
    tuple val(meta), path("*_expressionHist.png")          , emit: hist              , optional: true
    tuple val(meta), path("*_dupMatrix.txt")               , emit: dupmatrix         , optional: true
    tuple val(meta), path("*_intercept_slope.txt")         , emit: intercept_slope   , optional: true
    tuple val(meta), path("*_dup_intercept_mqc.txt")       , emit: multiqc_intercept , optional: true
    tuple val(meta), path("*_duprateExpDensCurve_mqc.txt") , emit: multiqc_curve     , optional: true
    // featureCounts / biotype outputs
    tuple val(meta), path("*.featureCounts.tsv")           , emit: featurecounts         , optional: true
    tuple val(meta), path("*.featureCounts.tsv.summary")   , emit: featurecounts_summary , optional: true
    tuple val(meta), path("*.biotype_counts.tsv")          , emit: biotype_counts_raw    , optional: true
    tuple val(meta), path("*.biotype_counts_mqc.tsv")      , emit: biotype_counts        , optional: true
    tuple val(meta), path("*.biotype_counts_rrna_mqc.tsv") , emit: biotype_rrna          , optional: true
    // RSeQC: bam_stat
    tuple val(meta), path("*.bam_stat.txt")                , emit: bamstat_txt                , optional: true
    // RSeQC: infer_experiment
    tuple val(meta), path("*.infer_experiment.txt")        , emit: inferexperiment_txt        , optional: true
    // RSeQC: read_duplication
    tuple val(meta), path("*.pos.DupRate.xls")             , emit: readduplication_pos_xls    , optional: true
    tuple val(meta), path("*.seq.DupRate.xls")             , emit: readduplication_seq_xls    , optional: true
    tuple val(meta), path("*.DupRate_plot.png")             , emit: readduplication_plot       , optional: true
    // RSeQC: read_distribution
    tuple val(meta), path("*.read_distribution.txt")       , emit: readdistribution_txt       , optional: true
    // RSeQC: junction_annotation
    tuple val(meta), path("*.junction.xls")                , emit: junctionannotation_xls     , optional: true
    tuple val(meta), path("*.junction.bed")                , emit: junctionannotation_bed     , optional: true
    tuple val(meta), path("*.junction_plot.r")             , emit: junctionannotation_rscript , optional: true
    tuple val(meta), path("*.junction_annotation.txt")     , emit: junctionannotation_log     , optional: true
    tuple val(meta), path("*.splice_events.png")           , emit: junctionannotation_events  , optional: true
    tuple val(meta), path("*.splice_junction.png")         , emit: junctionannotation_junctions, optional: true
    // RSeQC: junction_saturation
    tuple val(meta), path("*.junctionSaturation_plot.r")   , emit: junctionsaturation_rscript , optional: true
    tuple val(meta), path("*.junctionSaturation_plot.png") , emit: junctionsaturation_plot    , optional: true
    tuple val(meta), path("*.junctionSaturation_summary.txt"), emit: junctionsaturation_summary, optional: true
    // RSeQC: inner_distance
    tuple val(meta), path("*.inner_distance_freq.txt")     , emit: innerdistance_freq         , optional: true
    tuple val(meta), path("*.inner_distance.txt")          , emit: innerdistance_txt          , optional: true
    tuple val(meta), path("*.inner_distance_plot.png")     , emit: innerdistance_plot         , optional: true
    tuple val(meta), path("*.inner_distance_plot.r")       , emit: innerdistance_rscript      , optional: true
    tuple val(meta), path("*.inner_distance_summary.txt")  , emit: innerdistance_summary      , optional: true
    // versions
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def strandedness = 0
    if (meta.strandedness == 'forward') {
        strandedness = 1
    } else if (meta.strandedness == 'reverse') {
        strandedness = 2
    }
    def paired = meta.single_end ? '' : '--paired'
    """
    rustqc rna \\
        ${bam} \\
        --gtf ${gtf} \\
        --stranded ${strandedness} \\
        ${paired} \\
        --threads ${task.cpus} \\
        --outdir . \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rustqc: \$(rustqc --version 2>&1 | sed 's/rustqc //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_duprateExpDens.png
    touch ${prefix}_duprateExpBoxplot.png
    touch ${prefix}_expressionHist.png
    touch ${prefix}_dupMatrix.txt
    touch ${prefix}_intercept_slope.txt
    touch ${prefix}_dup_intercept_mqc.txt
    touch ${prefix}_duprateExpDensCurve_mqc.txt
    touch ${prefix}.featureCounts.tsv
    touch ${prefix}.featureCounts.tsv.summary
    touch ${prefix}.biotype_counts.tsv
    touch ${prefix}.biotype_counts_mqc.tsv
    touch ${prefix}.biotype_counts_rrna_mqc.tsv
    touch ${prefix}.bam_stat.txt
    touch ${prefix}.infer_experiment.txt
    touch ${prefix}.pos.DupRate.xls
    touch ${prefix}.seq.DupRate.xls
    touch ${prefix}.read_distribution.txt
    touch ${prefix}.junction.xls
    touch ${prefix}.junction.bed
    touch ${prefix}.junction_plot.r
    touch ${prefix}.junction_annotation.txt
    touch ${prefix}.splice_events.png
    touch ${prefix}.splice_junction.png
    touch ${prefix}.junctionSaturation_plot.r
    touch ${prefix}.junctionSaturation_plot.png
    touch ${prefix}.junctionSaturation_summary.txt
    touch ${prefix}.DupRate_plot.png
    touch ${prefix}.inner_distance.txt
    touch ${prefix}.inner_distance_freq.txt
    touch ${prefix}.inner_distance_plot.png
    touch ${prefix}.inner_distance_plot.r
    touch ${prefix}.inner_distance_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rustqc: \$(rustqc --version 2>&1 | sed 's/rustqc //')
    END_VERSIONS
    """
}
