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
    bam_bai       // channel: [ meta: Map, bam: Path, bai: Path ]
    bed           // file: genome.bed
    rseqc_modules // list: rseqc modules to run

    emit:
    bam_bai.map { meta, bam, bai ->
        //
        // Run RSeQC bam_stat.py
        //
        def bamstat = 'bam_stat' in rseqc_modules
            ? RSEQC_BAMSTAT(meta, bam)
            : null

        //
        // Run RSeQC inner_distance.py
        //
        def innerdistance = 'inner_distance' in rseqc_modules
            ? RSEQC_INNERDISTANCE(meta, bam, bed)
            : null

        //
        // Run RSeQC infer_experiment.py
        //
        def inferexperiment = 'infer_experiment' in rseqc_modules
            ? RSEQC_INFEREXPERIMENT(meta, bam, bed)
            : null

        //
        // Run RSeQC junction_annotation.py
        //
        def junctionannotation = 'junction_annotation' in rseqc_modules
            ? RSEQC_JUNCTIONANNOTATION(meta, bam, bed)
            : null

        //
        // Run RSeQC junction_saturation.py
        //
        def junctionsaturation = 'junction_saturation' in rseqc_modules
            ? RSEQC_JUNCTIONSATURATION(meta, bam, bed)
            : null

        //
        // Run RSeQC read_distribution.py
        //
        def readdistribution = 'read_distribution' in rseqc_modules
            ? RSEQC_READDISTRIBUTION(meta, bam, bed)
            : null

        //
        // Run RSeQC read_duplication.py
        //
        def readduplication = 'read_duplication' in rseqc_modules
            ? RSEQC_READDUPLICATION(meta, bam)
            : null

        //
        // Run RSeQC tin.py
        //
        def tin = 'tin' in rseqc_modules
            ? RSEQC_READDUPLICATION(meta, bam, bai, bed)
            : null

        [ meta, bamstat, innerdistance, inferexperiment, junctionannotation, junctionsaturation, readdistribution, readduplication, tin ]
    }

    // channel: [ meta: Map, txt ]
    // channel: [ meta: Map, {txt, pdf, r} ]
    // channel: [ meta: Map, txt ]
    // channel: [ meta: Map, {bed, xls, pdf, r, log} ]
    // channel: [ meta: Map, {pdf, r} ]
    // channel: [ meta: Map, txt ]
    // channel: [ meta: Map, {xls, pdf, r} ]
    // channel: [ meta: Map, txt ]

}
