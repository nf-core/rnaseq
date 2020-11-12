/*
 * Run RSeQC modules
 */

params.bamstat_options            = [:]
params.innerdistance_options      = [:]
params.inferexperiment_options    = [:]
params.junctionannotation_options = [:]
params.junctionsaturation_options = [:]
params.readdistribution_options   = [:]
params.readduplication_options    = [:]

include { RSEQC_BAMSTAT            } from '../software/rseqc/bamstat/main'            addParams( options: params.bamstat_options            )
include { RSEQC_INNERDISTANCE      } from '../software/rseqc/innerdistance/main'      addParams( options: params.innerdistance_options      )
include { RSEQC_INFEREXPERIMENT    } from '../software/rseqc/inferexperiment/main'    addParams( options: params.inferexperiment_options    )
include { RSEQC_JUNCTIONANNOTATION } from '../software/rseqc/junctionannotation/main' addParams( options: params.junctionannotation_options )
include { RSEQC_JUNCTIONSATURATION } from '../software/rseqc/junctionsaturation/main' addParams( options: params.junctionsaturation_options )
include { RSEQC_READDISTRIBUTION   } from '../software/rseqc/readdistribution/main'   addParams( options: params.readdistribution_options   )
include { RSEQC_READDUPLICATION    } from '../software/rseqc/readduplication/main'    addParams( options: params.readduplication_options    )

workflow RSEQC {
    take:
    bam           // channel: [ val(meta), [ ban ] ]
    bed           //    file: /path/to/genome.bed
    rseqc_modules //    list: rseqc modules to run

    main:
    /*
     * Run RSeQC bam_stat.py
     */
    version       = Channel.empty()
    bamstat_txt   = Channel.empty()
    if ('bam_stat' in rseqc_modules) {
        RSEQC_BAMSTAT ( bam )
        bamstat_txt = RSEQC_BAMSTAT.out.txt
        version     = RSEQC_BAMSTAT.out.version
    }

    /*
     * Run RSeQC inner_distance.py
     */
    innerdistance_distance = Channel.empty()
    innerdistance_freq     = Channel.empty()
    innerdistance_mean     = Channel.empty()
    innerdistance_pdf      = Channel.empty()
    innerdistance_rscript  = Channel.empty()
    if ('inner_distance' in rseqc_modules) {
        RSEQC_INNERDISTANCE ( bam, bed )
        innerdistance_distance = RSEQC_INNERDISTANCE.out.distance
        innerdistance_freq     = RSEQC_INNERDISTANCE.out.freq
        innerdistance_mean     = RSEQC_INNERDISTANCE.out.mean
        innerdistance_pdf      = RSEQC_INNERDISTANCE.out.pdf
        innerdistance_rscript  = RSEQC_INNERDISTANCE.out.rscript
        version                = RSEQC_INNERDISTANCE.out.version
    }

    /*
     * Run RSeQC infer_experiment.py
     */
    inferexperiment_txt   = Channel.empty()
    if ('infer_experiment' in rseqc_modules) {
        RSEQC_INFEREXPERIMENT ( bam, bed )
        inferexperiment_txt = RSEQC_INFEREXPERIMENT.out.txt
        version             = RSEQC_INFEREXPERIMENT.out.version
    }

    /*
     * Run RSeQC junction_annotation.py
     */
    junctionannotation_bed          = Channel.empty()
    junctionannotation_interact_bed = Channel.empty()
    junctionannotation_xls          = Channel.empty()
    junctionannotation_pdf          = Channel.empty()
    junctionannotation_events_pdf   = Channel.empty()
    junctionannotation_rscript      = Channel.empty()
    junctionannotation_log          = Channel.empty()
    if ('junction_annotation' in rseqc_modules) {
        RSEQC_JUNCTIONANNOTATION ( bam, bed )
        junctionannotation_bed          = RSEQC_JUNCTIONANNOTATION.out.bed
        junctionannotation_interact_bed = RSEQC_JUNCTIONANNOTATION.out.interact_bed
        junctionannotation_xls          = RSEQC_JUNCTIONANNOTATION.out.xls
        junctionannotation_pdf          = RSEQC_JUNCTIONANNOTATION.out.pdf
        junctionannotation_events_pdf   = RSEQC_JUNCTIONANNOTATION.out.events_pdf
        junctionannotation_rscript      = RSEQC_JUNCTIONANNOTATION.out.rscript
        junctionannotation_log          = RSEQC_JUNCTIONANNOTATION.out.log
        version                         = RSEQC_JUNCTIONANNOTATION.out.version
    }

    /*
     * Run RSeQC junction_saturation.py
     */
    junctionsaturation_pdf     = Channel.empty()
    junctionsaturation_rscript = Channel.empty()
    if ('junction_saturation' in rseqc_modules) {
        RSEQC_JUNCTIONSATURATION ( bam, bed )
        junctionsaturation_pdf     = RSEQC_JUNCTIONSATURATION.out.pdf
        junctionsaturation_rscript = RSEQC_JUNCTIONSATURATION.out.rscript
        version                    = RSEQC_JUNCTIONSATURATION.out.version
    }

    /*
     * Run RSeQC read_distribution.py
     */
    readdistribution_txt = Channel.empty()
    if ('read_distribution' in rseqc_modules) {
        RSEQC_READDISTRIBUTION ( bam, bed )
        readdistribution_txt = RSEQC_READDISTRIBUTION.out.txt
        version              = RSEQC_READDISTRIBUTION.out.version
    }

    /*
     * Run RSeQC read_duplication.py
     */
    readduplication_seq_xls = Channel.empty()
    readduplication_pos_xls = Channel.empty()
    readduplication_pdf     = Channel.empty()
    readduplication_rscript = Channel.empty()
    if ('read_duplication' in rseqc_modules) {
        RSEQC_READDUPLICATION ( bam )
        readduplication_seq_xls = RSEQC_READDUPLICATION.out.seq_xls
        readduplication_pos_xls = RSEQC_READDUPLICATION.out.pos_xls
        readduplication_pdf     = RSEQC_READDUPLICATION.out.pdf
        readduplication_rscript = RSEQC_READDUPLICATION.out.rscript
        version                 = RSEQC_READDUPLICATION.out.version
    }

    emit:
    bamstat_txt                     // channel: [ val(meta), txt ]

    innerdistance_distance          // channel: [ val(meta), txt ]
    innerdistance_freq              // channel: [ val(meta), txt ]
    innerdistance_mean              // channel: [ val(meta), txt ]
    innerdistance_pdf               // channel: [ val(meta), pdf ]
    innerdistance_rscript           // channel: [ val(meta), r   ]

    inferexperiment_txt             // channel: [ val(meta), txt ]

    junctionannotation_bed          // channel: [ val(meta), bed ]
    junctionannotation_interact_bed // channel: [ val(meta), bed ]
    junctionannotation_xls          // channel: [ val(meta), xls ]
    junctionannotation_pdf          // channel: [ val(meta), pdf ]
    junctionannotation_events_pdf   // channel: [ val(meta), pdf ]
    junctionannotation_rscript      // channel: [ val(meta), r   ]
    junctionannotation_log          // channel: [ val(meta), log ]

    junctionsaturation_pdf          // channel: [ val(meta), pdf ]
    junctionsaturation_rscript      // channel: [ val(meta), r   ]

    readdistribution_txt            // channel: [ val(meta), txt ]

    readduplication_seq_xls         // channel: [ val(meta), xls ]
    readduplication_pos_xls         // channel: [ val(meta), xls ]
    readduplication_pdf             // channel: [ val(meta), pdf ]
    readduplication_rscript         // channel: [ val(meta), r   ]

    version                         //    path: *.version.txt
}
