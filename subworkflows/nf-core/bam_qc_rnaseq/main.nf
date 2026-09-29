//
// Run post-alignment QC tools on RNA-seq BAM files
//

include { DUPRADAR                        } from '../../../modules/nf-core/dupradar/main'
include { PRESEQ_LCEXTRAP                 } from '../../../modules/nf-core/preseq/lcextrap/main'
include { QUALIMAP_RNASEQ                 } from '../../../modules/nf-core/qualimap/rnaseq/main'
include { SUBREAD_FEATURECOUNTS           } from '../../../modules/nf-core/subread/featurecounts/main'
include { CUSTOM_MULTIQCCUSTOMBIOTYPE     } from '../../../modules/nf-core/custom/multiqccustombiotype/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_QUALIMAP } from '../../../modules/nf-core/samtools/sort/main'
include { BAM_RSEQC                       } from '../bam_rseqc/main'
include { BamQcRnaseq                     } from './types'

workflow BAM_QC_RNASEQ {

    take:
    ch_bam_bai         // channel: [ val(meta), path(bam), path(bai) ]
    ch_gtf             // channel: [ val(meta), path(gtf) ]
    ch_gene_bed        // channel: path(bed)
    ch_fasta_fai       // channel: [ val(meta), path(fasta), path(fai) ]
    ch_biotypes_header // channel: [ val(meta), path(biotypes_header.txt) ]
    tools   // val(list)   - e.g. ['preseq', 'biotype_qc', 'qualimap', 'dupradar', 'rseqc_bam_stat', 'rseqc_infer_experiment', ...]
    biotype // val(string) - e.g. "gene_type" or "gene_biotype"

    main:
    def rseqc_modules = tools.findAll { tool -> tool.startsWith('rseqc_') }.collect { tool -> tool.replace('rseqc_', '') }

    ch_genome_bam = ch_bam_bai.map { meta, bam, _bai -> [ meta, bam ] }

    //
    // MODULE: Preseq library complexity
    //
    PRESEQ_LCEXTRAP (
        ch_genome_bam.filter { 'preseq' in tools }
    )

    //
    // MODULE: Feature biotype QC using featureCounts
    //
    SUBREAD_FEATURECOUNTS (
        ch_genome_bam
            .combine(ch_gtf.map { _meta, gtf -> gtf })
            .filter { 'biotype_qc' in tools && biotype }
    )

    ch_featurecounts_counts = SUBREAD_FEATURECOUNTS.out.map { r -> [r.meta, r.counts] }

    CUSTOM_MULTIQCCUSTOMBIOTYPE (
        ch_featurecounts_counts,
        ch_biotypes_header
    )
    ch_biotype_tsv  = CUSTOM_MULTIQCCUSTOMBIOTYPE.out.map { r -> [r.meta, r.tsv] }
    ch_biotype_rrna = CUSTOM_MULTIQCCUSTOMBIOTYPE.out.map { r -> [r.meta, r.rrna] }

    //
    // MODULE: Qualimap (name-sorted BAM via samtools sort)
    // Requires ext.args = '-n' to be set by the caller for SAMTOOLS_SORT_QUALIMAP
    //
    SAMTOOLS_SORT_QUALIMAP (
        ch_genome_bam.filter { 'qualimap' in tools },
        ch_fasta_fai,
        ''
    )

    QUALIMAP_RNASEQ (
        SAMTOOLS_SORT_QUALIMAP.out.bam,
        ch_gtf
    )

    //
    // MODULE: dupRadar
    //
    DUPRADAR (
        ch_genome_bam.filter { 'dupradar' in tools },
        ch_gtf
    )

    //
    // SUBWORKFLOW: RSeQC
    //
    BAM_RSEQC (
        ch_bam_bai
            .map { meta, bam, bai -> [ meta, [ bam, bai ] ] }
            .filter { rseqc_modules.size() > 0 },
        ch_gene_bed,
        rseqc_modules
    )

    ch_preseq_lc_extrap      = PRESEQ_LCEXTRAP.out.map { r -> [r.meta, r.lc_extrap] }
    ch_preseq_log            = PRESEQ_LCEXTRAP.out.map { r -> [r.meta, r.log] }
    ch_featurecounts_summary = SUBREAD_FEATURECOUNTS.out.map { r -> [r.meta, r.summary] }
    ch_dupradar_multiqc      = DUPRADAR.out.map { r -> [r.meta, r.multiqc] }

    // Each tool group that ran joins onto this per-sample skeleton; groups
    // that were skipped are left unset and become null, so no join is made
    // against a channel that is empty by construction.
    ch_results = ch_bam_bai.map { meta, _bam, _bai -> [meta.id, [:]] }

    if ('preseq' in tools) {
        // Remainder join: preseq can fail on low-duplication BAMs, and callers
        // may set errorStrategy 'ignore' rather than lose the sample.
        ch_results = ch_results
            .join(PRESEQ_LCEXTRAP.out.map { r -> [r.meta.id, r] }, by: [0], remainder: true)
            .map { id, fields, preseq -> [id, fields + [preseq: preseq]] }
    }

    if ('biotype_qc' in tools && biotype) {
        ch_results = ch_results
            .join(SUBREAD_FEATURECOUNTS.out.map { r -> [r.meta.id, r] }, by: [0])
            .join(CUSTOM_MULTIQCCUSTOMBIOTYPE.out.map { r -> [r.meta.id, record(tsv: r.tsv, rrna: r.rrna)] }, by: [0])
            .map { id, fields, featurecounts, biotype_files ->
                [id, fields + [
                    featurecounts: featurecounts,
                    biotype:       biotype_files
                ]]
            }
    }

    if ('qualimap' in tools) {
        ch_results = ch_results
            .join(QUALIMAP_RNASEQ.out.map { meta, dir -> [meta.id, dir] }, by: [0])
            .map { id, fields, dir -> [id, fields + [qualimap: dir]] }
    }

    if ('dupradar' in tools) {
        ch_results = ch_results
            .join(DUPRADAR.out.map { r -> [r.meta.id, r] }, by: [0])
            .map { id, fields, dupradar -> [id, fields + [dupradar: dupradar]] }
    }

    if (rseqc_modules.size() > 0) {
        ch_results = ch_results
            .join(BAM_RSEQC.out.results.map { r -> [r.id, r] }, by: [0])
            .map { id, fields, rseqc -> [id, fields + [rseqc: rseqc]] }
    }

    ch_results = ch_results.map { id, fields ->
        record(
            id:            id,
            preseq:        fields.preseq,
            featurecounts: fields.featurecounts,
            biotype:       fields.biotype,
            qualimap:      fields.qualimap,
            dupradar:      fields.dupradar,
            rseqc:         fields.rseqc
        )
    }

    // Aggregate MultiQC-compatible output files
    ch_multiqc_files = channel.empty()
        .mix(ch_preseq_lc_extrap)
        .mix(ch_biotype_tsv)
        .mix(QUALIMAP_RNASEQ.out)
        .mix(ch_dupradar_multiqc)
        .mix(BAM_RSEQC.out.bamstat_txt)
        .mix(BAM_RSEQC.out.inferexperiment_txt)
        .mix(BAM_RSEQC.out.innerdistance_freq)
        .mix(BAM_RSEQC.out.junctionannotation_log)
        .mix(BAM_RSEQC.out.junctionsaturation_rscript)
        .mix(BAM_RSEQC.out.readdistribution_txt)
        .mix(BAM_RSEQC.out.readduplication_pos_xls)
        .mix(BAM_RSEQC.out.tin_txt)

    // `remainder: true` needed because every contributor is gated behind
    // `tools` / `rseqc_modules`.
    ch_per_sample_mqc_bundle = ch_preseq_lc_extrap
        .join(ch_biotype_tsv,                              remainder: true)
        .join(QUALIMAP_RNASEQ.out,                         remainder: true)
        .join(ch_dupradar_multiqc,                         remainder: true)
        .join(BAM_RSEQC.out.bamstat_txt,                   remainder: true)
        .join(BAM_RSEQC.out.inferexperiment_txt,           remainder: true)
        .join(BAM_RSEQC.out.innerdistance_freq,            remainder: true)
        .join(BAM_RSEQC.out.junctionannotation_log,        remainder: true)
        .join(BAM_RSEQC.out.junctionsaturation_rscript,    remainder: true)
        .join(BAM_RSEQC.out.readdistribution_txt,          remainder: true)
        .join(BAM_RSEQC.out.readduplication_pos_xls,       remainder: true)
        .join(BAM_RSEQC.out.tin_txt,                       remainder: true)
        .map { row -> [row[0], row.drop(1).findAll { f -> f != null }.collectMany { e -> (e instanceof List) ? e : [e] }] }

    emit:
    // Aggregated
    multiqc_files = ch_multiqc_files // channel: [ val(meta), path(files) ]

    // Preseq
    preseq_lc_extrap = ch_preseq_lc_extrap // channel: [ val(meta), path(txt) ]
    preseq_log       = ch_preseq_log       // channel: [ val(meta), path(log) ]

    // Biotype QC
    featurecounts_counts  = ch_featurecounts_counts     // channel: [ val(meta), path(txt) ]
    featurecounts_summary = ch_featurecounts_summary    // channel: [ val(meta), path(txt) ]
    biotype_tsv           = ch_biotype_tsv              // channel: [ val(meta), path(tsv) ]
    biotype_rrna          = ch_biotype_rrna             // channel: [ val(meta), path(tsv) ]

    // Qualimap
    qualimap_results = QUALIMAP_RNASEQ.out // channel: [ val(meta), path(dir) ]

    // dupRadar
    dupradar_scatter2d       = DUPRADAR.out.map { r -> [r.meta, r.scatter2d] }       // channel: [ val(meta), path(pdf) ]
    dupradar_boxplot         = DUPRADAR.out.map { r -> [r.meta, r.boxplot] }         // channel: [ val(meta), path(pdf) ]
    dupradar_hist            = DUPRADAR.out.map { r -> [r.meta, r.hist] }            // channel: [ val(meta), path(pdf) ]
    dupradar_dupmatrix       = DUPRADAR.out.map { r -> [r.meta, r.dupmatrix] }       // channel: [ val(meta), path(txt) ]
    dupradar_intercept_slope = DUPRADAR.out.map { r -> [r.meta, r.intercept_slope] } // channel: [ val(meta), path(txt) ]
    dupradar_multiqc         = ch_dupradar_multiqc         // channel: [ val(meta), path(txt) ]

    // RSeQC
    inferexperiment_txt             = BAM_RSEQC.out.inferexperiment_txt             // channel: [ val(meta), path(txt) ]
    bamstat_txt                     = BAM_RSEQC.out.bamstat_txt                     // channel: [ val(meta), path(txt) ]
    innerdistance_all               = BAM_RSEQC.out.innerdistance_all               // channel: [ val(meta), path(txt/pdf/r) ]
    innerdistance_distance          = BAM_RSEQC.out.innerdistance_distance          // channel: [ val(meta), path(txt) ]
    innerdistance_freq              = BAM_RSEQC.out.innerdistance_freq              // channel: [ val(meta), path(txt) ]
    innerdistance_mean              = BAM_RSEQC.out.innerdistance_mean              // channel: [ val(meta), path(txt) ]
    innerdistance_pdf               = BAM_RSEQC.out.innerdistance_pdf               // channel: [ val(meta), path(pdf) ]
    innerdistance_rscript           = BAM_RSEQC.out.innerdistance_rscript           // channel: [ val(meta), path(r) ]
    junctionannotation_all          = BAM_RSEQC.out.junctionannotation_all          // channel: [ val(meta), path(bed/xls/pdf/r/log) ]
    junctionannotation_bed          = BAM_RSEQC.out.junctionannotation_bed          // channel: [ val(meta), path(bed) ]
    junctionannotation_interact_bed = BAM_RSEQC.out.junctionannotation_interact_bed // channel: [ val(meta), path(bed) ]
    junctionannotation_xls          = BAM_RSEQC.out.junctionannotation_xls          // channel: [ val(meta), path(xls) ]
    junctionannotation_pdf          = BAM_RSEQC.out.junctionannotation_pdf          // channel: [ val(meta), path(pdf) ]
    junctionannotation_events_pdf   = BAM_RSEQC.out.junctionannotation_events_pdf   // channel: [ val(meta), path(pdf) ]
    junctionannotation_rscript      = BAM_RSEQC.out.junctionannotation_rscript      // channel: [ val(meta), path(r) ]
    junctionannotation_log          = BAM_RSEQC.out.junctionannotation_log          // channel: [ val(meta), path(log) ]
    junctionsaturation_all          = BAM_RSEQC.out.junctionsaturation_all          // channel: [ val(meta), path(pdf/r) ]
    junctionsaturation_pdf          = BAM_RSEQC.out.junctionsaturation_pdf          // channel: [ val(meta), path(pdf) ]
    junctionsaturation_rscript      = BAM_RSEQC.out.junctionsaturation_rscript      // channel: [ val(meta), path(r) ]
    readdistribution_txt            = BAM_RSEQC.out.readdistribution_txt            // channel: [ val(meta), path(txt) ]
    readduplication_all             = BAM_RSEQC.out.readduplication_all             // channel: [ val(meta), path(xls/pdf/r) ]
    readduplication_seq_xls         = BAM_RSEQC.out.readduplication_seq_xls         // channel: [ val(meta), path(xls) ]
    readduplication_pos_xls         = BAM_RSEQC.out.readduplication_pos_xls         // channel: [ val(meta), path(xls) ]
    readduplication_pdf             = BAM_RSEQC.out.readduplication_pdf             // channel: [ val(meta), path(pdf) ]
    readduplication_rscript         = BAM_RSEQC.out.readduplication_rscript         // channel: [ val(meta), path(r) ]
    tin_txt                         = BAM_RSEQC.out.tin_txt                         // channel: [ val(meta), path(txt) ]
    per_sample_mqc_bundle           = ch_per_sample_mqc_bundle                      // channel: [ val(meta), list(files) ]
    results                         = ch_results                                    // channel: BamQcRnaseq

}
