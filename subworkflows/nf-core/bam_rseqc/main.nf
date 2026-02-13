nextflow.preview.types = true

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

include { InnerDistanceResult      } from '../../../modules/nf-core/rseqc/innerdistance/main'
include { JunctionAnnotationResult } from '../../../modules/nf-core/rseqc/junctionannotation/main'
include { ReadDuplicationResult    } from '../../../modules/nf-core/rseqc/readduplication/main'

record RSeQCResult {
    meta:                    Map
    bamstat:                 Path?
    inferexperiment:         Path?
    junction_annotation:     JunctionAnnotationResult?
    junctionsaturation_pdf:  Path?
    junctionsaturation_r:    Path?
    read_duplication:        ReadDuplicationResult?
    readdistribution:        Path?
    inner_distance:          InnerDistanceResult?
    tin:                     Path?
}

workflow BAM_RSEQC {
    take:
    bam_bai       // channel: [ val(meta), [ bam, bai ] ]
    bed           // channel: [ genome.bed ]
    rseqc_modules //    list: rseqc modules to run

    main:

    bam = bam_bai.map{ it -> [ it[0], it[1][0], it[1][1] ] }

    // Per-module channels - records for multi-output modules, paths for single-output
    bamstat_txt            = channel.empty()
    inferexperiment_txt    = channel.empty()
    ch_junction_annotation = channel.empty()  // [meta, JunctionAnnotationResult]
    junctionsaturation_pdf = channel.empty()
    junctionsaturation_r   = channel.empty()
    ch_read_duplication    = channel.empty()  // [meta, ReadDuplicationResult]
    readdistribution_txt   = channel.empty()
    ch_inner_distance      = channel.empty()  // [meta, InnerDistanceResult]
    tin_txt                = channel.empty()

    if ('bam_stat' in rseqc_modules) {
        RSEQC_BAMSTAT(bam)
        bamstat_txt = RSEQC_BAMSTAT.out.map { r -> [r.meta, r.txt] }
    }

    if ('inner_distance' in rseqc_modules) {
        RSEQC_INNERDISTANCE(bam, bed)
        ch_inner_distance = RSEQC_INNERDISTANCE.out.map { r -> [r.meta, r] }
    }

    if ('infer_experiment' in rseqc_modules) {
        RSEQC_INFEREXPERIMENT(bam, bed)
        inferexperiment_txt = RSEQC_INFEREXPERIMENT.out.map { r -> [r.meta, r.txt] }
    }

    if ('junction_annotation' in rseqc_modules) {
        RSEQC_JUNCTIONANNOTATION(bam, bed)
        ch_junction_annotation = RSEQC_JUNCTIONANNOTATION.out.map { r -> [r.meta, r] }
    }

    if ('junction_saturation' in rseqc_modules) {
        RSEQC_JUNCTIONSATURATION(bam, bed)
        junctionsaturation_pdf = RSEQC_JUNCTIONSATURATION.out.map { r -> [r.meta, r.pdf] }
        junctionsaturation_r   = RSEQC_JUNCTIONSATURATION.out.map { r -> [r.meta, r.rscript] }
    }

    if ('read_distribution' in rseqc_modules) {
        RSEQC_READDISTRIBUTION(bam, bed)
        readdistribution_txt = RSEQC_READDISTRIBUTION.out.map { r -> [r.meta, r.txt] }
    }

    if ('read_duplication' in rseqc_modules) {
        RSEQC_READDUPLICATION(bam)
        ch_read_duplication = RSEQC_READDUPLICATION.out.map { r -> [r.meta, r] }
    }

    if ('tin' in rseqc_modules) {
        RSEQC_TIN(bam, bed)
        tin_txt = RSEQC_TIN.out.map { r -> [r.meta, r.txt] }
    }

    // Join all per-module outputs by meta, using bam as the driver channel
    // (ensures all samples are represented even when modules are skipped).
    // Modules that weren't run produce null fields via remainder: true.
    emit:
    result = bam
        .map { meta, _b, _i -> [meta] }
        .join(bamstat_txt,            by: [0], remainder: true)
        .join(inferexperiment_txt,    by: [0], remainder: true)
        .join(ch_junction_annotation, by: [0], remainder: true)
        .join(junctionsaturation_pdf, by: [0], remainder: true)
        .join(junctionsaturation_r,   by: [0], remainder: true)
        .join(ch_read_duplication,    by: [0], remainder: true)
        .join(readdistribution_txt,   by: [0], remainder: true)
        .join(ch_inner_distance,      by: [0], remainder: true)
        .join(tin_txt,                by: [0], remainder: true)
        .map { meta, bamstat, inferexp, junction_annot, juncsat_pdf, juncsat_r, read_dup, readdist, inner_dist, tin ->
            record(
                meta: meta,
                bamstat: bamstat, inferexperiment: inferexp,
                junction_annotation: junction_annot,
                junctionsaturation_pdf: juncsat_pdf, junctionsaturation_r: juncsat_r,
                read_duplication: read_dup,
                readdistribution: readdist,
                inner_distance: inner_dist,
                tin: tin
            )
        }
    inferexperiment = inferexperiment_txt
    multiqc_files = bamstat_txt
        .mix(inferexperiment_txt)
        .mix(ch_inner_distance.map { meta, r -> [meta, r.freq] })
        .mix(ch_junction_annotation.map { meta, r -> [meta, r.log] })
        .mix(junctionsaturation_r)
        .mix(readdistribution_txt)
        .mix(ch_read_duplication.map { meta, r -> [meta, r.pos_xls] })
        .mix(tin_txt)
        .map { _meta, file -> file }
}
