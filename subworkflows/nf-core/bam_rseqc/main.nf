nextflow.enable.types = true

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
include { Bam; Rseqc } from '../../../modules/nf-core/types'

workflow BAM_RSEQC {
    take:
    ch_bam_bai: Channel<Bam>
    bed: Value<Path>
    rseqc_modules: List<String>

    main:
    // Every field starts null and is overwritten by the join of the tool that ran, so
    // skipped tools leave a null field and no join is made against an empty channel.
    ch_results = ch_bam_bai.map { r ->
        record(
            id:                 r.id,
            meta:               r.meta,
            bamstat:            null,
            inferexperiment:    null,
            innerdistance:      null,
            junctionannotation: null,
            junctionsaturation: null,
            readdistribution:   null,
            readduplication:    null,
            tin:                null
        )
    }

    if ('bam_stat' in rseqc_modules) {
        ch_results = ch_results.join(RSEQC_BAMSTAT(ch_bam_bai), by: 'id')
    }

    if ('inner_distance' in rseqc_modules) {
        ch_results = ch_results.join(
            RSEQC_INNERDISTANCE(ch_bam_bai, bed).map { r -> record(id: r.id, innerdistance: r) },
            by: 'id'
        )
    }

    if ('infer_experiment' in rseqc_modules) {
        ch_results = ch_results.join(RSEQC_INFEREXPERIMENT(ch_bam_bai, bed), by: 'id')
    }

    if ('junction_annotation' in rseqc_modules) {
        ch_results = ch_results.join(
            RSEQC_JUNCTIONANNOTATION(ch_bam_bai, bed).map { r -> record(id: r.id, junctionannotation: r) },
            by: 'id'
        )
    }

    if ('junction_saturation' in rseqc_modules) {
        ch_results = ch_results.join(
            RSEQC_JUNCTIONSATURATION(ch_bam_bai, bed).map { r -> record(id: r.id, junctionsaturation: r) },
            by: 'id'
        )
    }

    if ('read_distribution' in rseqc_modules) {
        ch_results = ch_results.join(RSEQC_READDISTRIBUTION(ch_bam_bai, bed), by: 'id')
    }

    if ('read_duplication' in rseqc_modules) {
        ch_results = ch_results.join(
            RSEQC_READDUPLICATION(ch_bam_bai).map { r -> record(id: r.id, readduplication: r) },
            by: 'id'
        )
    }

    if ('tin' in rseqc_modules) {
        ch_results = ch_results.join(
            RSEQC_TIN(ch_bam_bai, bed).map { r -> record(id: r.id, tin: r) },
            by: 'id'
        )
    }

    emit:
    ch_results
}
