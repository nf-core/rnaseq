nextflow.enable.types = true

//
// umicollapse, index BAM file and run samtools stats, flagstat and idxstats
//

include { UMICOLLAPSE        } from '../../../modules/nf-core/umicollapse/main'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS } from '../bam_stats_samtools/main'
include { BamInput; UmicollapseDedupBam } from '../../../modules/nf-core/types'
include { SamtoolsIndexResult } from '../../../modules/nf-core/samtools/index/main'
include { UmicollapseResult } from '../../../modules/nf-core/umicollapse/main'

//
// UMICollapse options
//
def umicollapseArgs(grouping_method: String?, umi_separator: String?, single_end: Boolean) -> String {
    def algos = ['directional': 'dir', 'adjacency': 'adj', 'cluster': 'cc']
    def paired_args = single_end ? '' : '--paired --remove-unpaired --remove-chimeric'
    def algo_arg = grouping_method ? "--algo '${algos[grouping_method] ?: ''}'" : ''
    def separator_arg = umi_separator ? "--umi-sep '${umi_separator}'" : ''
    return ['--two-pass', paired_args, algo_arg, separator_arg].join(' ').trim()
}

record UmicollapseArgs {
    samtools_index: String?
}

workflow BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE {
    take:
    ch_bam_bai: Channel<BamInput>
    tool_args: UmicollapseArgs // samtools index options
    grouping_method: String? // UMI grouping method
    umi_separator: String? // UMI separator

    main:
    //
    // umicollapse in bam mode (thus hardcode mode input to 'bam')
    //
    def ch_dedup: Channel<UmicollapseResult> = UMICOLLAPSE(ch_bam_bai.map { r -> r + record(args: umicollapseArgs(grouping_method, umi_separator, r.meta.single_end)) }, 'bam')
        .filter { r -> r.bam != null }

    //
    // Index BAM file and run samtools stats, flagstat and idxstats
    //
    def ch_index: Channel<SamtoolsIndexResult> = SAMTOOLS_INDEX(ch_dedup.map { r -> r + record(args: tool_args.samtools_index) })
    ch_indexed = ch_dedup.join(ch_index, by: 'id')

    ch_results = ch_indexed.join(BAM_STATS_SAMTOOLS(ch_indexed, null, null), by: 'id')

    emit:
    ch_results
}
