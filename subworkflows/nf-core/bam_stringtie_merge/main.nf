nextflow.enable.types = true

include { STRINGTIE_STRINGTIE } from '../../../modules/nf-core/stringtie/stringtie/main'
include { STRINGTIE_MERGE     } from '../../../modules/nf-core/stringtie/merge/main'
include { StringtieInput; StringtieMerged } from './types'

workflow BAM_STRINGTIE_MERGE {
    take:
    ch_bams: Channel<StringtieInput>
    mode: Value<List<String>>
    chrgtf: Value<Tuple<Map, Path?>>

    main:

    ch_assemblies = STRINGTIE_STRINGTIE(
        ch_bams,
        mode,
        chrgtf.map { _meta, gtf -> gtf }
    )

    ch_to_merge = ch_assemblies
        .collect()
        .map { assemblies ->
            record(
                id:         'stringtie_merge',
                meta:       [id: 'stringtie_merge'],
                gtf:        assemblies.collect { r -> r.transcript_gtf }.toSorted { f -> f.name },
                assemblies: assemblies
            )
        }

    ch_merged = STRINGTIE_MERGE(ch_to_merge, chrgtf)

    ch_results = ch_to_merge.combine(merged_gtf: ch_merged.map { r -> r.merged_gtf })

    emit:
    ch_results
}
