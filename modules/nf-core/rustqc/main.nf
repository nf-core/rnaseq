process RUSTQC {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2a/2a8a0514855c54307399fd0f664c2685e76c8cc07631e767c1e37c575b18d59f/data'
        : 'community.wave.seqera.io/library/rustqc:0.2.1--00df1502b490e005'}"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(gtf)

    output:
    tuple val(meta), path("${prefix}/dupradar/*"),                                                      emit: dupradar
    tuple val(meta), path("${prefix}/featurecounts/*"),                                                 emit: featurecounts
    tuple val(meta), path("${prefix}/preseq/*"),                                                        emit: preseq
    tuple val(meta), path("${prefix}/samtools/*"),                                                      emit: samtools
    tuple val(meta), path("${prefix}/rseqc/**"),                                                        emit: rseqc
    tuple val(meta), path("${prefix}/qualimap/**"),                                                     emit: qualimap
    tuple val("${task.process}"), val('rustqc'), eval("rustqc --version 2>&1 | sed -n '1s/rustqc //; 1s/ .*//p'"),  emit: versions_rustqc, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def paired = meta.single_end ? '' : '--paired'
    """
    rustqc rna \\
        ${bam} \\
        --gtf ${gtf} \\
        ${paired} \\
        --threads ${task.cpus} \\
        --outdir ${prefix} \\
        --sample-name ${prefix} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/{dupradar,featurecounts,preseq,samtools} \\
            ${prefix}/rseqc/{bam_stat,infer_experiment,read_duplication,read_distribution,junction_annotation,junction_saturation,inner_distance,tin} \\
            ${prefix}/qualimap/{raw_data_qualimapReport,images_qualimapReport}

    touch ${prefix}/dupradar/${prefix}_{duprateExpDens,duprateExpBoxplot,expressionHist}.png \\
          ${prefix}/dupradar/${prefix}_{dupMatrix,intercept_slope,dup_intercept_mqc,duprateExpDensCurve_mqc}.txt

    touch ${prefix}/featurecounts/${prefix}.featureCounts.tsv \\
          ${prefix}/featurecounts/${prefix}.featureCounts.tsv.summary \\
          ${prefix}/featurecounts/${prefix}.featureCounts.biotype.tsv.summary \\
          ${prefix}/featurecounts/${prefix}.{biotype_counts,biotype_counts_mqc,biotype_counts_rrna_mqc}.tsv

    touch ${prefix}/rseqc/bam_stat/${prefix}.bam_stat.txt \\
          ${prefix}/rseqc/infer_experiment/${prefix}.infer_experiment.txt \\
          ${prefix}/rseqc/read_duplication/${prefix}.{pos.DupRate,seq.DupRate}.xls \\
          ${prefix}/rseqc/read_duplication/${prefix}.DupRate_plot.{r,png} \\
          ${prefix}/rseqc/read_distribution/${prefix}.read_distribution.txt \\
          ${prefix}/rseqc/junction_annotation/${prefix}.{junction.xls,junction.bed,junction_plot.r,junction_annotation.log,splice_events.png,splice_junction.png} \\
          ${prefix}/rseqc/junction_saturation/${prefix}.junctionSaturation_{plot.r,plot.png,summary.txt} \\
          ${prefix}/rseqc/inner_distance/${prefix}.inner_distance{.txt,_freq.txt,_plot.r,_plot.png,_summary.txt,_mean.txt} \\
          ${prefix}/rseqc/tin/${prefix}.{tin.xls,summary.txt}

    touch ${prefix}/preseq/${prefix}.lc_extrap.txt \\
          ${prefix}/samtools/${prefix}.{flagstat,idxstats,stats} \\
          ${prefix}/qualimap/rnaseq_qc_results.txt \\
          ${prefix}/qualimap/qualimapReport.html
    touch "${prefix}/qualimap/raw_data_qualimapReport/coverage_profile_along_genes_(total).txt"
    """
}
