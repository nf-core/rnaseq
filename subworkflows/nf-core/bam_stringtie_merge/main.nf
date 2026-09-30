nextflow.enable.types = true

include { STRINGTIE_STRINGTIE } from '../../../modules/nf-core/stringtie/stringtie/main'
include { STRINGTIE_MERGE     } from '../../../modules/nf-core/stringtie/merge/main'
include { StringtieInput; StringtieMerged } from '../../../modules/nf-core/types'
include { StringtieMergeResult } from '../../../modules/nf-core/stringtie/merge/main'
include { StringtieResult } from '../../../modules/nf-core/stringtie/stringtie/main'

workflow BAM_STRINGTIE_MERGE {
    take:
    ch_bams: Channel<StringtieInput>
    mode: Value<List<String>>
    chrgtf: Value<Path?>

    main:

    def ch_assemblies: Channel<StringtieResult> = STRINGTIE_STRINGTIE(
        ch_bams,
        mode,
        chrgtf
    )

    ch_to_merge = ch_assemblies
        .collect()
        .flatMap { assemblies ->
            assemblies.isEmpty() ? [] : [
                record(
                    id:         'stringtie_merge',
                    meta:       [id: 'stringtie_merge'],
                    gtf:        assemblies.collect { r -> r.transcript_gtf }.toSorted { f -> f.name },
                    assemblies: assemblies
                )
            ]
        }

    def ch_merged: Channel<StringtieMergeResult> = STRINGTIE_MERGE(ch_to_merge, chrgtf)

    ch_results = ch_to_merge.join(ch_merged.map { r -> record(id: r.id, merged_gtf: r.merged_gtf) }, by: 'id')

    emit:
    ch_results
}
