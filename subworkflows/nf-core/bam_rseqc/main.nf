//
// Run RSeQC modules
//

include { RSEQC_BAMSTAT            } from '../../../modules/nf-core/rseqc/bamstat/main'
include { RSEQC_INNERDISTANCE      } from '../../../modules/nf-core/rseqc/innerdistance/main'
include { RSEQC_INFEREXPERIMENT    } from '../../../modules/nf-core/rseqc/inferexperiment/main'
include { RSEQC_JUNCTIONANNOTATION } from '../../../modules/nf-core/rseqc/junctionannotation/main'
include { RSEQC_JUNCTIONSATURATION } from '../../../modules/nf-core/rseqc/junctionsaturation/main'
include { RSEQC_READDISTRIBUTION   } from '../../../modules/nf-core/rseqc/readdistribution/main'
include { RSEQC_READDUPLICATION    } from '../../../modules/nf-core/rseqc/readduplication/main'
include { RSEQC_TIN                } from '../../../modules/nf-core/rseqc/tin/main'
include { Rseqc                    } from './types'

workflow BAM_RSEQC {
    take:
    bam_bai       // channel: [ val(meta), [ bam, bai ] ]
    bed           // channel: [ genome.bed ]
    rseqc_modules //    list: rseqc modules to run

    main:

    // bam = bam_bai.map{ [ it[0], it[1][0], it[1][1] ] }
    bam = bam_bai.map{ it -> [ it[0], it[1][0], it[1][1] ] }

    // Each tool that runs joins its outputs onto this per-sample skeleton;
    // fields for tools that did not run are left unset and become null in
    // the final record, so no join is ever made against an empty channel.
    ch_results = bam.map { meta, _bam, _bai -> [meta.id, [:]] }

    //
    // Run RSeQC bam_stat.py
    //
    bamstat_txt = channel.empty()

    if ('bam_stat' in rseqc_modules) {
        RSEQC_BAMSTAT(bam)
        bamstat_txt = RSEQC_BAMSTAT.out.txt
        ch_results = ch_results
            .join(bamstat_txt.map { meta, txt -> [meta.id, txt] }, by: [0])
            .map { id, fields, txt -> [id, fields + [bamstat: txt]] }
    }

    //
    // Run RSeQC inner_distance.py
    //
    innerdistance_all      = channel.empty()
    innerdistance_distance = channel.empty()
    innerdistance_freq     = channel.empty()
    innerdistance_mean     = channel.empty()
    innerdistance_pdf      = channel.empty()
    innerdistance_rscript  = channel.empty()

    if ('inner_distance' in rseqc_modules) {
        RSEQC_INNERDISTANCE(bam, bed)
        innerdistance_distance = RSEQC_INNERDISTANCE.out.distance
        innerdistance_freq     = RSEQC_INNERDISTANCE.out.freq
        innerdistance_mean     = RSEQC_INNERDISTANCE.out.mean
        innerdistance_pdf      = RSEQC_INNERDISTANCE.out.pdf
        innerdistance_rscript  = RSEQC_INNERDISTANCE.out.rscript
        innerdistance_all      = innerdistance_distance.mix(innerdistance_freq, innerdistance_mean, innerdistance_pdf, innerdistance_rscript)

        // Single-end samples only produce the distance file
        ch_results = ch_results
            .join(innerdistance_distance.map { meta, f -> [meta.id, f] }, by: [0])
            .join(innerdistance_freq.map { meta, f -> [meta.id, f] }, by: [0], remainder: true)
            .join(innerdistance_mean.map { meta, f -> [meta.id, f] }, by: [0], remainder: true)
            .join(innerdistance_pdf.map { meta, f -> [meta.id, f] }, by: [0], remainder: true)
            .join(innerdistance_rscript.map { meta, f -> [meta.id, f] }, by: [0], remainder: true)
            .map { id, fields, distance, freq, mean, pdf, rscript ->
                [id, fields + [innerdistance: record(distance: distance, freq: freq, mean: mean, pdf: pdf, rscript: rscript)]]
            }
    }

    //
    // Run RSeQC infer_experiment.py
    //
    inferexperiment_txt = channel.empty()
    if ('infer_experiment' in rseqc_modules) {
        RSEQC_INFEREXPERIMENT(bam, bed)
        inferexperiment_txt = RSEQC_INFEREXPERIMENT.out.txt
        ch_results = ch_results
            .join(inferexperiment_txt.map { meta, txt -> [meta.id, txt] }, by: [0])
            .map { id, fields, txt -> [id, fields + [inferexperiment: txt]] }
    }

    //
    // Run RSeQC junction_annotation.py
    //
    junctionannotation_all          = channel.empty()
    junctionannotation_bed          = channel.empty()
    junctionannotation_interact_bed = channel.empty()
    junctionannotation_xls          = channel.empty()
    junctionannotation_pdf          = channel.empty()
    junctionannotation_events_pdf   = channel.empty()
    junctionannotation_rscript      = channel.empty()
    junctionannotation_log          = channel.empty()

    if ('junction_annotation' in rseqc_modules) {
        RSEQC_JUNCTIONANNOTATION(bam, bed)
        junctionannotation_bed          = RSEQC_JUNCTIONANNOTATION.out.bed
        junctionannotation_interact_bed = RSEQC_JUNCTIONANNOTATION.out.interact_bed
        junctionannotation_xls          = RSEQC_JUNCTIONANNOTATION.out.xls
        junctionannotation_pdf          = RSEQC_JUNCTIONANNOTATION.out.pdf
        junctionannotation_events_pdf   = RSEQC_JUNCTIONANNOTATION.out.events_pdf
        junctionannotation_rscript      = RSEQC_JUNCTIONANNOTATION.out.rscript
        junctionannotation_log          = RSEQC_JUNCTIONANNOTATION.out.log
        junctionannotation_all          = junctionannotation_bed.mix(junctionannotation_interact_bed, junctionannotation_xls, junctionannotation_pdf, junctionannotation_events_pdf, junctionannotation_rscript, junctionannotation_log)

        ch_results = ch_results
            .join(junctionannotation_xls.map { meta, f -> [meta.id, f] }, by: [0])
            .join(junctionannotation_rscript.map { meta, f -> [meta.id, f] }, by: [0])
            .join(junctionannotation_log.map { meta, f -> [meta.id, f] }, by: [0])
            .join(junctionannotation_bed.map { meta, f -> [meta.id, f] }, by: [0], remainder: true)
            .join(junctionannotation_interact_bed.map { meta, f -> [meta.id, f] }, by: [0], remainder: true)
            .join(junctionannotation_pdf.map { meta, f -> [meta.id, f] }, by: [0], remainder: true)
            .join(junctionannotation_events_pdf.map { meta, f -> [meta.id, f] }, by: [0], remainder: true)
            .map { id, fields, xls, rscript, log, bed_file, interact_bed, pdf, events_pdf ->
                [id, fields + [junctionannotation: record(bed: bed_file, interact_bed: interact_bed, xls: xls, pdf: pdf, events_pdf: events_pdf, rscript: rscript, log: log)]]
            }
    }

    //
    // Run RSeQC junction_saturation.py
    //
    junctionsaturation_all     = channel.empty()
    junctionsaturation_pdf     = channel.empty()
    junctionsaturation_rscript = channel.empty()

    if ('junction_saturation' in rseqc_modules) {
        RSEQC_JUNCTIONSATURATION(bam, bed)
        junctionsaturation_pdf     = RSEQC_JUNCTIONSATURATION.out.pdf
        junctionsaturation_rscript = RSEQC_JUNCTIONSATURATION.out.rscript
        junctionsaturation_all     = junctionsaturation_pdf.mix(junctionsaturation_rscript)

        ch_results = ch_results
            .join(junctionsaturation_pdf.map { meta, f -> [meta.id, f] }, by: [0])
            .join(junctionsaturation_rscript.map { meta, f -> [meta.id, f] }, by: [0])
            .map { id, fields, pdf, rscript ->
                [id, fields + [junctionsaturation: record(pdf: pdf, rscript: rscript)]]
            }
    }

    //
    // Run RSeQC read_distribution.py
    //
    readdistribution_txt = channel.empty()

    if ('read_distribution' in rseqc_modules) {
        RSEQC_READDISTRIBUTION(bam, bed)
        readdistribution_txt = RSEQC_READDISTRIBUTION.out.txt
        ch_results = ch_results
            .join(readdistribution_txt.map { meta, txt -> [meta.id, txt] }, by: [0])
            .map { id, fields, txt -> [id, fields + [readdistribution: txt]] }
    }

    //
    // Run RSeQC read_duplication.py
    //
    readduplication_all     = channel.empty()
    readduplication_seq_xls = channel.empty()
    readduplication_pos_xls = channel.empty()
    readduplication_pdf     = channel.empty()
    readduplication_rscript = channel.empty()

    if ('read_duplication' in rseqc_modules) {
        RSEQC_READDUPLICATION(bam )
        readduplication_seq_xls = RSEQC_READDUPLICATION.out.seq_xls
        readduplication_pos_xls = RSEQC_READDUPLICATION.out.pos_xls
        readduplication_pdf     = RSEQC_READDUPLICATION.out.pdf
        readduplication_rscript = RSEQC_READDUPLICATION.out.rscript
        readduplication_all     = readduplication_seq_xls.mix(readduplication_pos_xls, readduplication_pdf, readduplication_rscript)

        ch_results = ch_results
            .join(readduplication_seq_xls.map { meta, f -> [meta.id, f] }, by: [0])
            .join(readduplication_pos_xls.map { meta, f -> [meta.id, f] }, by: [0])
            .join(readduplication_pdf.map { meta, f -> [meta.id, f] }, by: [0])
            .join(readduplication_rscript.map { meta, f -> [meta.id, f] }, by: [0])
            .map { id, fields, seq_xls, pos_xls, pdf, rscript ->
                [id, fields + [readduplication: record(seq_xls: seq_xls, pos_xls: pos_xls, pdf: pdf, rscript: rscript)]]
            }
    }

    //
    // Run RSeQC tin.py
    //
    tin_txt = channel.empty()

    if ('tin' in rseqc_modules) {
        RSEQC_TIN(bam, bed)
        tin_txt      = RSEQC_TIN.out.txt
        ch_results = ch_results
            .join(tin_txt.map { meta, txt -> [meta.id, txt] }, by: [0])
            .join(RSEQC_TIN.out.xls.map { meta, xls -> [meta.id, xls] }, by: [0])
            .map { id, fields, txt, xls -> [id, fields + [tin: record(txt: txt, xls: xls)]] }
    }

    ch_results = ch_results.map { id, fields ->
        record(
            id:                 id,
            bamstat:            fields.bamstat,
            inferexperiment:    fields.inferexperiment,
            innerdistance:      fields.innerdistance,
            junctionannotation: fields.junctionannotation,
            junctionsaturation: fields.junctionsaturation,
            readdistribution:   fields.readdistribution,
            readduplication:    fields.readduplication,
            tin:                fields.tin
        )
    }

    emit:
    bamstat_txt                     // channel: [ val(meta), txt ]

    innerdistance_all               // channel: [ val(meta), {txt, pdf, r} ]
    innerdistance_distance          // channel: [ val(meta), txt ]
    innerdistance_freq              // channel: [ val(meta), txt ]
    innerdistance_mean              // channel: [ val(meta), txt ]
    innerdistance_pdf               // channel: [ val(meta), pdf ]
    innerdistance_rscript           // channel: [ val(meta), r   ]

    inferexperiment_txt             // channel: [ val(meta), txt ]

    junctionannotation_all          // channel: [ val(meta), {bed, xls, pdf, r, log} ]
    junctionannotation_bed          // channel: [ val(meta), bed ]
    junctionannotation_interact_bed // channel: [ val(meta), bed ]
    junctionannotation_xls          // channel: [ val(meta), xls ]
    junctionannotation_pdf          // channel: [ val(meta), pdf ]
    junctionannotation_events_pdf   // channel: [ val(meta), pdf ]
    junctionannotation_rscript      // channel: [ val(meta), r   ]
    junctionannotation_log          // channel: [ val(meta), log ]

    junctionsaturation_all          // channel: [ val(meta), {pdf, r} ]
    junctionsaturation_pdf          // channel: [ val(meta), pdf ]
    junctionsaturation_rscript      // channel: [ val(meta), r   ]

    readdistribution_txt            // channel: [ val(meta), txt ]

    readduplication_all             // channel: [ val(meta), {xls, pdf, r} ]
    readduplication_seq_xls         // channel: [ val(meta), xls ]
    readduplication_pos_xls         // channel: [ val(meta), xls ]
    readduplication_pdf             // channel: [ val(meta), pdf ]
    readduplication_rscript         // channel: [ val(meta), r   ]

    tin_txt                         // channel: [ val(meta), txt ]

    results = ch_results            // channel: Rseqc
}
