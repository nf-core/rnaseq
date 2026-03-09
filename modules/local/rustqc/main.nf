process RUSTQC {
    tag "$meta.id"
    label 'process_high'

    container "ghcr.io/seqeralabs/rustqc:dev"

    input:
    tuple val(meta), path(bam), path(bai)
    path gtf

    output:
    tuple val(meta), path("rustqc")                                                , emit: results
    tuple val(meta), path("rustqc/rseqc/infer_experiment/*.infer_experiment.txt")  , emit: inferexperiment_txt, optional: true
    path "versions.yml"                                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def strandedness = 0
    if (meta.strandedness == 'forward') {
        strandedness = 1
    } else if (meta.strandedness == 'reverse') {
        strandedness = 2
    }
    def paired       = meta.single_end ? '' : '--paired'
    def biotype_attr = params.gencode ? "--biotype-attribute gene_type" : ''
    """
    rustqc rna \\
        ${bam} \\
        --gtf ${gtf} \\
        --stranded ${strandedness} \\
        ${paired} \\
        --threads ${task.cpus} \\
        --outdir rustqc \\
        ${biotype_attr} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rustqc: \$(rustqc --version | sed 's/rustqc //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p rustqc/dupradar
    mkdir -p rustqc/featurecounts
    mkdir -p rustqc/rseqc/bam_stat
    mkdir -p rustqc/rseqc/infer_experiment
    mkdir -p rustqc/rseqc/read_duplication
    mkdir -p rustqc/rseqc/read_distribution
    mkdir -p rustqc/rseqc/junction_annotation
    mkdir -p rustqc/rseqc/junction_saturation
    mkdir -p rustqc/rseqc/inner_distance
    mkdir -p rustqc/rseqc/tin
    mkdir -p rustqc/preseq
    mkdir -p rustqc/samtools
    mkdir -p rustqc/qualimap/raw_data_qualimapReport
    mkdir -p rustqc/qualimap/images_qualimapReport

    touch rustqc/dupradar/${prefix}_duprateExpDens.png
    touch rustqc/dupradar/${prefix}_duprateExpBoxplot.png
    touch rustqc/dupradar/${prefix}_expressionHist.png
    touch rustqc/dupradar/${prefix}_dupMatrix.txt
    touch rustqc/dupradar/${prefix}_intercept_slope.txt
    touch rustqc/dupradar/${prefix}_dup_intercept_mqc.txt
    touch rustqc/dupradar/${prefix}_duprateExpDensCurve_mqc.txt

    touch rustqc/featurecounts/${prefix}.featureCounts.tsv
    touch rustqc/featurecounts/${prefix}.featureCounts.tsv.summary
    touch rustqc/featurecounts/${prefix}.biotype_counts.tsv
    touch rustqc/featurecounts/${prefix}.biotype_counts_mqc.tsv
    touch rustqc/featurecounts/${prefix}.biotype_counts_rrna_mqc.tsv

    touch rustqc/rseqc/bam_stat/${prefix}.bam_stat.txt
    touch rustqc/rseqc/infer_experiment/${prefix}.infer_experiment.txt
    touch rustqc/rseqc/read_duplication/${prefix}.pos.DupRate.xls
    touch rustqc/rseqc/read_duplication/${prefix}.seq.DupRate.xls
    touch rustqc/rseqc/read_duplication/${prefix}.DupRate_plot.r
    touch rustqc/rseqc/read_duplication/${prefix}.DupRate_plot.png
    touch rustqc/rseqc/read_distribution/${prefix}.read_distribution.txt
    touch rustqc/rseqc/junction_annotation/${prefix}.junction.xls
    touch rustqc/rseqc/junction_annotation/${prefix}.junction.bed
    touch rustqc/rseqc/junction_annotation/${prefix}.junction_plot.r
    touch rustqc/rseqc/junction_annotation/${prefix}.junction_annotation.txt
    touch rustqc/rseqc/junction_annotation/${prefix}.splice_events.png
    touch rustqc/rseqc/junction_annotation/${prefix}.splice_junction.png
    touch rustqc/rseqc/junction_saturation/${prefix}.junctionSaturation_plot.r
    touch rustqc/rseqc/junction_saturation/${prefix}.junctionSaturation_plot.png
    touch rustqc/rseqc/junction_saturation/${prefix}.junctionSaturation_summary.txt
    touch rustqc/rseqc/inner_distance/${prefix}.inner_distance.txt
    touch rustqc/rseqc/inner_distance/${prefix}.inner_distance_freq.txt
    touch rustqc/rseqc/inner_distance/${prefix}.inner_distance_plot.r
    touch rustqc/rseqc/inner_distance/${prefix}.inner_distance_plot.png
    touch rustqc/rseqc/inner_distance/${prefix}.inner_distance_summary.txt
    touch rustqc/rseqc/inner_distance/${prefix}.inner_distance_mean.txt
    touch rustqc/rseqc/tin/${prefix}.tin.xls
    touch rustqc/rseqc/tin/${prefix}.summary.txt

    touch rustqc/preseq/${prefix}.lc_extrap.txt
    touch rustqc/samtools/${prefix}.flagstat
    touch rustqc/samtools/${prefix}.idxstats
    touch rustqc/samtools/${prefix}.stats

    touch rustqc/qualimap/rnaseq_qc_results.txt
    touch "rustqc/qualimap/raw_data_qualimapReport/coverage_profile_along_genes_(total).txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rustqc: \$(rustqc --version | sed 's/rustqc //')
    END_VERSIONS
    """
}
