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
    ch_bam_bai    // channel: [ val(meta), [ bam, bai ] ]
    ch_bed        //    file: /path/to/genome.bed
    rseqc_modules //    list: rseqc modules to run

    main:

    ch_versions = Channel.empty()

    ch_bam_bai
        .map { [ it[0], it[1] ] }
        .set { ch_bam }

    //
    // Run RSeQC bam_stat.py
    //
    ch_bamstat = Channel.empty()
    if ('bam_stat' in rseqc_modules) {
        RSEQC_BAMSTAT ( ch_bam )
        ch_bamstat  = RSEQC_BAMSTAT.out.txt
        ch_versions = ch_versions.mix(RSEQC_BAMSTAT.out.versions.first())
    }

    //
    // Run RSeQC inner_distance.py
    //
    ch_innerdistance_distance = Channel.empty()
    ch_innerdistance_freq     = Channel.empty()
    ch_innerdistance_mean     = Channel.empty()
    ch_innerdistance_pdf      = Channel.empty()
    ch_innerdistance_rscript  = Channel.empty()
    if ('inner_distance' in rseqc_modules) {
        RSEQC_INNERDISTANCE ( ch_bam, ch_bed )
        ch_innerdistance_distance = RSEQC_INNERDISTANCE.out.distance
        ch_innerdistance_freq     = RSEQC_INNERDISTANCE.out.freq
        ch_innerdistance_mean     = RSEQC_INNERDISTANCE.out.mean
        ch_innerdistance_pdf      = RSEQC_INNERDISTANCE.out.pdf
        ch_innerdistance_rscript  = RSEQC_INNERDISTANCE.out.rscript
        ch_inner_distance         = ch_innerdistance_distance.mix(ch_innerdistance_freq, ch_innerdistance_mean, ch_innerdistance_pdf, ch_innerdistance_rscript)
        ch_versions               = ch_versions.mix(RSEQC_INNERDISTANCE.out.versions.first())
    }

    //
    // Run RSeQC infer_experiment.py
    //
    ch_inferexperiment = Channel.empty()
    if ('infer_experiment' in rseqc_modules) {
        RSEQC_INFEREXPERIMENT ( ch_bam, ch_bed )
        ch_inferexperiment = RSEQC_INFEREXPERIMENT.out.txt
        ch_versions        = ch_versions.mix(RSEQC_INFEREXPERIMENT.out.versions.first())
    }

    //
    // Run RSeQC junction_annotation.py
    //
    ch_junctionannotation_bed          = Channel.empty()
    ch_junctionannotation_interact_bed = Channel.empty()
    ch_junctionannotation_xls          = Channel.empty()
    ch_junctionannotation_pdf          = Channel.empty()
    ch_junctionannotation_events_pdf   = Channel.empty()
    ch_junctionannotation_rscript      = Channel.empty()
    ch_junctionannotation_log          = Channel.empty()
    if ('junction_annotation' in rseqc_modules) {
        RSEQC_JUNCTIONANNOTATION ( ch_bam, ch_bed )
        ch_junctionannotation_bed          = RSEQC_JUNCTIONANNOTATION.out.bed
        ch_junctionannotation_interact_bed = RSEQC_JUNCTIONANNOTATION.out.interact_bed
        ch_junctionannotation_xls          = RSEQC_JUNCTIONANNOTATION.out.xls
        ch_junctionannotation_pdf          = RSEQC_JUNCTIONANNOTATION.out.pdf
        ch_junctionannotation_events_pdf   = RSEQC_JUNCTIONANNOTATION.out.events_pdf
        ch_junctionannotation_rscript      = RSEQC_JUNCTIONANNOTATION.out.rscript
        ch_junctionannotation_log          = RSEQC_JUNCTIONANNOTATION.out.log
        ch_junction_annotation             = ch_junctionannotation_bed.mix(ch_junctionannotation_interact_bed, ch_junctionannotation_xls, ch_junctionannotation_pdf, ch_junctionannotation_events_pdf, ch_junctionannotation_rscript, ch_junctionannotation_log)
        ch_versions                        = ch_versions.mix(RSEQC_JUNCTIONANNOTATION.out.versions.first())
    }

    //
    // Run RSeQC junction_saturation.py
    //
    ch_junctionsaturation_pdf     = Channel.empty()
    ch_junctionsaturation_rscript = Channel.empty()
    if ('junction_saturation' in rseqc_modules) {
        RSEQC_JUNCTIONSATURATION ( ch_bam, ch_bed )
        ch_junctionsaturation_pdf     = RSEQC_JUNCTIONSATURATION.out.pdf
        ch_junctionsaturation_rscript = RSEQC_JUNCTIONSATURATION.out.rscript
        ch_junction_saturation        = ch_junctionsaturation_pdf.mix(ch_junctionsaturation_rscript)
        ch_versions                   = ch_versions.mix(RSEQC_JUNCTIONSATURATION.out.versions.first())
    }

    //
    // Run RSeQC read_distribution.py
    //
    ch_readdistribution = Channel.empty()
    if ('read_distribution' in rseqc_modules) {
        RSEQC_READDISTRIBUTION ( ch_bam, ch_bed )
        ch_readdistribution = RSEQC_READDISTRIBUTION.out.txt
        ch_versions         = ch_versions.mix(RSEQC_READDISTRIBUTION.out.versions.first())
    }

    //
    // Run RSeQC read_duplication.py
    //
    ch_readduplication_seq_xls = Channel.empty()
    ch_readduplication_pos_xls = Channel.empty()
    ch_readduplication_pdf     = Channel.empty()
    ch_readduplication_rscript = Channel.empty()
    if ('read_duplication' in rseqc_modules) {
        RSEQC_READDUPLICATION ( ch_bam )
        ch_readduplication_seq_xls = RSEQC_READDUPLICATION.out.seq_xls
        ch_readduplication_pos_xls = RSEQC_READDUPLICATION.out.pos_xls
        ch_readduplication_pdf     = RSEQC_READDUPLICATION.out.pdf
        ch_readduplication_rscript = RSEQC_READDUPLICATION.out.rscript
        ch_read_duplication        = ch_readduplication_seq_xls.mix(ch_readduplication_pos_xls, ch_readduplication_pdf, ch_readduplication_rscript)
        ch_versions                = ch_versions.mix(RSEQC_READDUPLICATION.out.versions.first())
    }

    //
    // Run RSeQC tin.py
    //
    ch_tin = Channel.empty()
    if ('tin' in rseqc_modules) {
        RSEQC_TIN ( ch_bam_bai, ch_bed )
        ch_tin      = RSEQC_TIN.out.txt
        ch_versions = ch_versions.mix(RSEQC_TIN.out.versions.first())
    }

    emit:
    ch_bamstat                     // channel: [ val(meta), txt ]

    ch_inner_distance
    ch_innerdistance_distance          // channel: [ val(meta), txt ]
    ch_innerdistance_freq              // channel: [ val(meta), txt ]
    ch_innerdistance_mean              // channel: [ val(meta), txt ]
    ch_innerdistance_pdf               // channel: [ val(meta), pdf ]
    ch_innerdistance_rscript           // channel: [ val(meta), r   ]

    ch_inferexperiment             // channel: [ val(meta), txt ]

    ch_junction_annotation
    ch_junctionannotation_bed          // channel: [ val(meta), bed ]
    ch_junctionannotation_interact_bed // channel: [ val(meta), bed ]
    ch_junctionannotation_xls          // channel: [ val(meta), xls ]
    ch_junctionannotation_pdf          // channel: [ val(meta), pdf ]
    ch_junctionannotation_events_pdf   // channel: [ val(meta), pdf ]
    ch_junctionannotation_rscript      // channel: [ val(meta), r   ]
    ch_junctionannotation_log          // channel: [ val(meta), log ]

    ch_junction_saturation
    ch_junctionsaturation_pdf          // channel: [ val(meta), pdf ]
    ch_junctionsaturation_rscript      // channel: [ val(meta), r   ]

    ch_readdistribution            // channel: [ val(meta), txt ]

    ch_read_duplication
    ch_readduplication_seq_xls         // channel: [ val(meta), xls ]
    ch_readduplication_pos_xls         // channel: [ val(meta), xls ]
    ch_readduplication_pdf             // channel: [ val(meta), pdf ]
    ch_readduplication_rscript         // channel: [ val(meta), r   ]

    ch_tin                         // channel: [ val(meta), txt ]

    versions = ch_versions          // channel: [ versions.yml ]
}
