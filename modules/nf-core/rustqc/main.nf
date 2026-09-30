nextflow.enable.types = true

include { BamBaiInput } from '../types'

record RustqcResult {
    id:            String
    meta:          Map
    samtools:      RustqcSamtools
    preseq:        RustqcPreseq
    dupradar:      RustqcDupradar
    featurecounts: RustqcFeaturecounts
    biotype:       RustqcBiotype
    rseqc:         RustqcRseqc
    qualimap:      Path?
    all_files:     Set<Path>
}

record RustqcSamtools {
    stats:    Path?
    flagstat: Path?
    idxstats: Path?
}

record RustqcPreseq {
    lc_extrap: Path?
}

record RustqcDupradar {
    scatter2d:       Set<Path>
    boxplot:         Set<Path>
    hist:            Set<Path>
    dupmatrix:       Path?
    intercept_slope: Path?
    multiqc:         Set<Path>
}

record RustqcFeaturecounts {
    counts:  Path?
    summary: Path?
}

record RustqcBiotype {
    tsv:  Path?
    mqc:  Path?
    rrna: Path?
}

record RustqcTin {
    txt: Path?
    xls: Path?
}

record RustqcInnerdistance {
    distance: Path?
    freq:     Path?
    mean:     Path?
    summary:  Path?
    plot:     Set<Path>
    rscript:  Path?
}

record RustqcJunctionannotation {
    bed:          Path?
    interact_bed: Path?
    xls:          Path?
    log:          Path?
    plot:         Set<Path>
    rscript:      Path?
}

record RustqcJunctionsaturation {
    summary: Path?
    plot:    Set<Path>
    rscript: Path?
}

record RustqcReadduplication {
    seq_xls: Path?
    pos_xls: Path?
    plot:    Set<Path>
    rscript: Path?
}

record RustqcRseqc {
    bamstat:            Path?
    inferexperiment:    Path?
    readdistribution:   Path?
    tin:                RustqcTin
    innerdistance:      RustqcInnerdistance
    junctionannotation: RustqcJunctionannotation
    junctionsaturation: RustqcJunctionsaturation
    readduplication:    RustqcReadduplication
}

process RUSTQC {
    tag "$sample.meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2a/2a8a0514855c54307399fd0f664c2685e76c8cc07631e767c1e37c575b18d59f/data'
        : 'community.wave.seqera.io/library/rustqc:0.2.1--00df1502b490e005'}"

    input:
    sample: BamBaiInput
    gtf: Path

    output:
    record(
        id: sample.id,
        meta: sample.meta,
        samtools: record(
            stats:    file("*.stats", optional: true),
            flagstat: file("*.flagstat", optional: true),
            idxstats: file("*.idxstats", optional: true)
        ),
        preseq: record(
            lc_extrap: file("*.lc_extrap.txt", optional: true)
        ),
        dupradar: record(
            scatter2d:       files("*duprateExpDens.*", optional: true),
            boxplot:         files("*duprateExpBoxplot.*", optional: true),
            hist:            files("*expressionHist.*", optional: true),
            dupmatrix:       file("*dupMatrix.*", optional: true),
            intercept_slope: file("*intercept_slope.*", optional: true),
            multiqc:         files("*{dup_intercept,duprateExpDensCurve}_mqc.txt", optional: true)
        ),
        featurecounts: record(
            counts:  file("*.featureCounts.tsv", optional: true),
            summary: file("*.featureCounts.biotype.tsv.summary", optional: true)
        ),
        biotype: record(
            tsv:  file("*.biotype_counts.tsv", optional: true),
            mqc:  file("*.biotype_counts_mqc.tsv", optional: true),
            rrna: file("*.biotype_counts_rrna_mqc.tsv", optional: true)
        ),
        rseqc: record(
            bamstat:          file("*.bam_stat.txt", optional: true),
            inferexperiment:  file("*.infer_experiment.txt", optional: true),
            readdistribution: file("*.read_distribution.txt", optional: true),
            tin: record(
                txt: file("*.summary.txt", optional: true),
                xls: file("*.tin.xls", optional: true)
            ),
            innerdistance: record(
                distance: file("*.inner_distance.txt", optional: true),
                freq:     file("*.inner_distance_freq.txt", optional: true),
                mean:     file("*.inner_distance_mean.txt", optional: true),
                summary:  file("*.inner_distance_summary.txt", optional: true),
                plot:     files("*.inner_distance_plot.{png,svg}", optional: true),
                rscript:  file("*.inner_distance_plot.r", optional: true)
            ),
            junctionannotation: record(
                bed:          file("*.junction.bed", optional: true),
                interact_bed: file("*.junction.Interact.bed", optional: true),
                xls:          file("*.junction.xls", optional: true),
                log:          file("*.junction_annotation.log", optional: true),
                plot:         files("*.splice_*.{png,svg}", optional: true),
                rscript:      file("*.junction_plot.r", optional: true)
            ),
            junctionsaturation: record(
                summary: file("*.junctionSaturation_summary.txt", optional: true),
                plot:    files("*.junctionSaturation_plot.{png,svg}", optional: true),
                rscript: file("*.junctionSaturation_plot.r", optional: true)
            ),
            readduplication: record(
                seq_xls: file("*.seq.DupRate.xls", optional: true),
                pos_xls: file("*.pos.DupRate.xls", optional: true),
                plot:    files("*.DupRate_plot.{png,svg}", optional: true),
                rscript: file("*.DupRate_plot.r", optional: true)
            )
        ),
        qualimap:  file("${prefix}", optional: true),
        all_files: files("{*.txt,*.tsv,*.xls,*.log,*.stats,*.flagstat,*.idxstats,*.html,*_mqc.*,${prefix}/**}", optional: true)
    ) as RustqcResult

    topic:
    tuple(task.process, 'rustqc', eval("rustqc --version 2>&1 | sed -n '1s/rustqc //; 1s/ .*//p'")) >> 'versions'

    script:
    prefix = task.ext.prefix ?: "${sample.meta.id}"
    def args = task.ext.args ?: ''
    def paired = sample.meta.single_end ? '' : '--paired'
    // Flatten the tool subdirectories into the task root so every record field is a plain file,
    // keeping qualimap's own directory structure under the sample-named directory (the tool's
    // --outdir stays ${prefix} because featureCounts records it in its output header).
    def flatten = """
    dups=\$(find ${prefix} -type f -not -path '${prefix}/qualimap/*' | sed 's#.*/##' | sort | uniq -d)
    if [ -n "\$dups" ]; then echo "RustQC output basenames collide: \$dups" >&2; exit 1; fi
    find ${prefix} -type f -not -path '${prefix}/qualimap/*' -exec mv {} . \\;
    mv ${prefix}/qualimap qualimap_tmp
    rm -rf ${prefix}
    mv qualimap_tmp ${prefix}
    """
    """
    rustqc rna \\
        ${sample.bam} \\
        --gtf ${gtf} \\
        ${paired} \\
        --threads ${task.cpus} \\
        --outdir ${prefix} \\
        --sample-name ${prefix} \\
        ${args}
    ${flatten}
    """

    stub:
    prefix = task.ext.prefix ?: "${sample.meta.id}"
    // Flatten the tool subdirectories into the task root so every record field is a plain file,
    // keeping qualimap's own directory structure under the sample-named directory (the tool's
    // --outdir stays ${prefix} because featureCounts records it in its output header).
    def flatten = """
    dups=\$(find ${prefix} -type f -not -path '${prefix}/qualimap/*' | sed 's#.*/##' | sort | uniq -d)
    if [ -n "\$dups" ]; then echo "RustQC output basenames collide: \$dups" >&2; exit 1; fi
    find ${prefix} -type f -not -path '${prefix}/qualimap/*' -exec mv {} . \\;
    mv ${prefix}/qualimap qualimap_tmp
    rm -rf ${prefix}
    mv qualimap_tmp ${prefix}
    """
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
    ${flatten}
    """
}
