process RUSTQC {
    tag "$meta.id"
    label 'process_medium'

    container "ghcr.io/ewels/rustqc:latest"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(gtf)

    output:
    // dupRadar outputs
    tuple val(meta), path("*_duprateExpDens.png")          , emit: scatter2d
    tuple val(meta), path("*_duprateExpBoxplot.png")       , emit: boxplot
    tuple val(meta), path("*_expressionHist.png")          , emit: hist
    tuple val(meta), path("*_dupMatrix.txt")               , emit: dupmatrix
    tuple val(meta), path("*_intercept_slope.txt")         , emit: intercept_slope
    tuple val(meta), path("*_dup_intercept_mqc.txt")       , emit: multiqc_intercept
    tuple val(meta), path("*_duprateExpDensCurve_mqc.txt") , emit: multiqc_curve
    // featureCounts / biotype outputs
    tuple val(meta), path("*.featureCounts.tsv")           , emit: featurecounts
    tuple val(meta), path("*.featureCounts.tsv.summary")   , emit: featurecounts_summary
    tuple val(meta), path("*.biotype_counts.tsv")          , emit: biotype_counts_raw
    tuple val(meta), path("*.biotype_counts_mqc.tsv")      , emit: biotype_counts
    tuple val(meta), path("*.biotype_counts_rrna_mqc.tsv") , emit: biotype_rrna
    // MultiQC and versions
    tuple val(meta), path("*_mqc.{txt,tsv}")               , emit: multiqc
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
        ${gtf} \\
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rustqc: \$(rustqc --version 2>&1 | sed 's/rustqc //')
    END_VERSIONS
    """
}
