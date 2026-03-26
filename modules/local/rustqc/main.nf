process RUSTQC {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    // TODO: pin to a release tag before merge
    container "ghcr.io/seqeralabs/rustqc:dev"

    input:
    tuple val(meta), path(bam), path(bai)
    path gtf

    output:
    tuple val(meta), path("rustqc")                                                , emit: results
    tuple val(meta), path("rustqc/rseqc/infer_experiment/*.infer_experiment.txt")  , emit: inferexperiment_txt, optional: true
    tuple val("${task.process}"), val('rustqc'), eval("rustqc --version | sed 's/rustqc //'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def strandedness = meta.strandedness ?: 'unstranded'
    def paired       = meta.single_end ? '' : '--paired'
    """
    rustqc rna \\
        ${bam} \\
        --gtf ${gtf} \\
        --stranded ${strandedness} \\
        ${paired} \\
        --threads ${task.cpus} \\
        --outdir rustqc \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p rustqc/{dupradar,featurecounts,preseq,samtools} \
            rustqc/rseqc/{bam_stat,infer_experiment,read_duplication,read_distribution,junction_annotation,junction_saturation,inner_distance,tin} \
            rustqc/qualimap/{raw_data_qualimapReport,images_qualimapReport}

    touch rustqc/dupradar/${prefix}_{duprateExpDens,duprateExpBoxplot,expressionHist}.png \
          rustqc/dupradar/${prefix}_{dupMatrix,intercept_slope,dup_intercept_mqc,duprateExpDensCurve_mqc}.txt

    touch rustqc/featurecounts/${prefix}.featureCounts.tsv \
          rustqc/featurecounts/${prefix}.featureCounts.tsv.summary \
          rustqc/featurecounts/${prefix}.{biotype_counts,biotype_counts_mqc,biotype_counts_rrna_mqc}.tsv

    touch rustqc/rseqc/bam_stat/${prefix}.bam_stat.txt \
          rustqc/rseqc/infer_experiment/${prefix}.infer_experiment.txt \
          rustqc/rseqc/read_duplication/${prefix}.{pos.DupRate,seq.DupRate}.xls \
          rustqc/rseqc/read_duplication/${prefix}.DupRate_plot.{r,png} \
          rustqc/rseqc/read_distribution/${prefix}.read_distribution.txt \
          rustqc/rseqc/junction_annotation/${prefix}.{junction.xls,junction.bed,junction_plot.r,junction_annotation.txt,splice_events.png,splice_junction.png} \
          rustqc/rseqc/junction_saturation/${prefix}.junctionSaturation_{plot.r,plot.png,summary.txt} \
          rustqc/rseqc/inner_distance/${prefix}.inner_distance{.txt,_freq.txt,_plot.r,_plot.png,_summary.txt,_mean.txt} \
          rustqc/rseqc/tin/${prefix}.{tin.xls,summary.txt}

    touch rustqc/preseq/${prefix}.lc_extrap.txt \
          rustqc/samtools/${prefix}.{flagstat,idxstats,stats} \
          rustqc/qualimap/rnaseq_qc_results.txt
    touch "rustqc/qualimap/raw_data_qualimapReport/coverage_profile_along_genes_(total).txt"
    """
}
