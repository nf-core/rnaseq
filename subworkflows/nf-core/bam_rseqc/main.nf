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

workflow BAM_RSEQC {
    take:
    bam_bai       // channel: [ val(meta), [ bam, bai ] ]
    bed           // channel: [ genome.bed ]
    rseqc_modules //    list: rseqc modules to run

    main:

    bam = bam_bai.map{ [ it[0], it[1] ] }

    versions = Channel.empty()

    //
    // Run RSeQC bam_stat.py
    //
    bamstat_txt = Channel.empty()

    if ('bam_stat' in rseqc_modules) {
        RSEQC_BAMSTAT(bam)
        bamstat_txt = RSEQC_BAMSTAT.out.txt
        versions    = versions.mix(RSEQC_BAMSTAT.out.versions.first())
    }

    //
    // Run RSeQC inner_distance.py
    //
    innerdistance_all      = Channel.empty()
    innerdistance_distance = Channel.empty()
    innerdistance_freq     = Channel.empty()
    innerdistance_mean     = Channel.empty()
    innerdistance_pdf      = Channel.empty()
    innerdistance_rscript  = Channel.empty()

    if ('inner_distance' in rseqc_modules) {
        RSEQC_INNERDISTANCE(bam, bed)
        innerdistance_distance = RSEQC_INNERDISTANCE.out.distance
        innerdistance_freq     = RSEQC_INNERDISTANCE.out.freq
        innerdistance_mean     = RSEQC_INNERDISTANCE.out.mean
        innerdistance_pdf      = RSEQC_INNERDISTANCE.out.pdf
        innerdistance_rscript  = RSEQC_INNERDISTANCE.out.rscript
        innerdistance_all      = innerdistance_distance.mix(innerdistance_freq, innerdistance_mean, innerdistance_pdf, innerdistance_rscript)
        versions               = versions.mix(RSEQC_INNERDISTANCE.out.versions.first())
    }

    //
    // Run RSeQC infer_experiment.py
    //
    inferexperiment_txt = Channel.empty()
    if ('infer_experiment' in rseqc_modules) {
        RSEQC_INFEREXPERIMENT(bam, bed)
        inferexperiment_txt = RSEQC_INFEREXPERIMENT.out.txt
        versions            = versions.mix(RSEQC_INFEREXPERIMENT.out.versions.first())
    }

    //
    // Run RSeQC junction_annotation.py
    //
    junctionannotation_all          = Channel.empty()
    junctionannotation_bed          = Channel.empty()
    junctionannotation_interact_bed = Channel.empty()
    junctionannotation_xls          = Channel.empty()
    junctionannotation_pdf          = Channel.empty()
    junctionannotation_events_pdf   = Channel.empty()
    junctionannotation_rscript      = Channel.empty()
    junctionannotation_log          = Channel.empty()

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
        versions                        = versions.mix(RSEQC_JUNCTIONANNOTATION.out.versions.first())
    }

    //
    // Run RSeQC junction_saturation.py
    //
    junctionsaturation_all     = Channel.empty()
    junctionsaturation_pdf     = Channel.empty()
    junctionsaturation_rscript = Channel.empty()

    if ('junction_saturation' in rseqc_modules) {
        RSEQC_JUNCTIONSATURATION(bam, bed)
        junctionsaturation_pdf     = RSEQC_JUNCTIONSATURATION.out.pdf
        junctionsaturation_rscript = RSEQC_JUNCTIONSATURATION.out.rscript
        junctionsaturation_all     = junctionsaturation_pdf.mix(junctionsaturation_rscript)
        versions                   = versions.mix(RSEQC_JUNCTIONSATURATION.out.versions.first())
    }

    //
    // Run RSeQC read_distribution.py
    //
    readdistribution_txt = Channel.empty()

    if ('read_distribution' in rseqc_modules) {
        RSEQC_READDISTRIBUTION(bam, bed)
        readdistribution_txt = RSEQC_READDISTRIBUTION.out.txt
        versions            = versions.mix(RSEQC_READDISTRIBUTION.out.versions.first())
    }

    //
    // Run RSeQC read_duplication.py
    //
    readduplication_all     = Channel.empty()
    readduplication_seq_xls = Channel.empty()
    readduplication_pos_xls = Channel.empty()
    readduplication_pdf     = Channel.empty()
    readduplication_rscript = Channel.empty()

    if ('read_duplication' in rseqc_modules) {
        RSEQC_READDUPLICATION(bam )
        readduplication_seq_xls = RSEQC_READDUPLICATION.out.seq_xls
        readduplication_pos_xls = RSEQC_READDUPLICATION.out.pos_xls
        readduplication_pdf     = RSEQC_READDUPLICATION.out.pdf
        readduplication_rscript = RSEQC_READDUPLICATION.out.rscript
        readduplication_all     = readduplication_seq_xls.mix(readduplication_pos_xls, readduplication_pdf, readduplication_rscript)
        versions                = versions.mix(RSEQC_READDUPLICATION.out.versions.first())
    }

    //
    // Run RSeQC tin.py
    //
    tin_txt = Channel.empty()

    if ('tin' in rseqc_modules) {
        RSEQC_TIN(bam_bai, bed)
        tin_txt      = RSEQC_TIN.out.txt
        versions    = versions.mix(RSEQC_TIN.out.versions.first())
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

    versions                        // channel: [ versions.yml ]
}
