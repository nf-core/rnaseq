//
// BAM deduplication with UMI processing
//

include { BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE as BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE_TRANSCRIPTOME } from '../bam_dedup_stats_samtools_umicollapse'
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME       } from '../bam_dedup_stats_samtools_umitools'
include { BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE as BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE_GENOME        } from '../bam_dedup_stats_samtools_umicollapse'
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME              } from '../bam_dedup_stats_samtools_umitools'
include { BAM_SORT_STATS_SAMTOOLS                                                                    } from '../bam_sort_stats_samtools'

include { UMITOOLS_PREPAREFORRSEM                                                                    } from '../../../modules/nf-core/umitools/prepareforrsem'
include { SAMTOOLS_SORT                                                                              } from '../../../modules/nf-core/samtools/sort/main'

record SamtoolsStatsFiles {
    stats:    Path
    flagstat: Path
    idxstats: Path
}

record UmitoolsDedupStats {
    edit_distance:    Path
    per_umi:          Path
    umi_per_position: Path
}

record UmiDedupTranscriptome {
    bam:              Path
    dedup_bam:        Path
    sorted_bam:       Path
    sorted_bam_index: Path
    filtered_bam:     Path?
    stats:            Path
    flagstat:         Path
    idxstats:         Path
    tsv:              UmitoolsDedupStats?
}

record UmiDedupBam {
    id:                       String
    meta:                     Map
    bam:                      Path
    bai:                      Path
    genomic_dedup_log:        Path
    transcriptomic_dedup_log: Path?
    prepare_for_rsem_log:     Path?
    genome:                   SamtoolsStatsFiles
    transcriptome:            UmiDedupTranscriptome?
    tsv:                      UmitoolsDedupStats?
}

workflow BAM_DEDUP_UMI {
    take:
    ch_genome_bam // channel: [ val(meta), path(bam), path(bai) ]
    ch_fasta_fai // channel: [ val(meta), path(fasta), path(fai) ]
    umi_dedup_tool // string: 'umicollapse' or 'umitools'
    umitools_dedup_stats // boolean: whether to generate UMI-tools dedup stats
    ch_transcriptome_bam // channel: [ val(meta), path(bam) ]
    ch_transcript_fasta_fai // channel: [ val(meta), path(fasta), path(fai) ]
    umitools_dedup_primary_only // boolean: whether to filter to primary alignments before dedup

    main:
    ch_tsv_edit_distance = channel.empty()
    ch_tsv_per_umi = channel.empty()
    ch_tsv_umi_per_position = channel.empty()
    ch_genomic_dedup_log = channel.empty()
    ch_transcriptomic_dedup_log = channel.empty()
    ch_genome_dedup = channel.empty()
    ch_transcriptome_dedup = channel.empty()

    if (umi_dedup_tool != "umicollapse" && umi_dedup_tool != "umitools") {
        error("Unknown umi_dedup_tool '${umi_dedup_tool}'")
    }

    // Genome BAM deduplication
    if (umi_dedup_tool == "umicollapse") {
        BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE_GENOME(
            ch_genome_bam
        )
        UMI_DEDUP_GENOME = BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE_GENOME
        ch_genomic_dedup_log = UMI_DEDUP_GENOME.out.dedup_stats
        ch_genome_dedup = UMI_DEDUP_GENOME.out.results.map { r -> [r.id, r.bam, r.bai, r.dedup_stats, r.samtools, null] }
    }
    else if (umi_dedup_tool == "umitools") {
        BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME(
            ch_genome_bam,
            umitools_dedup_stats,
            umitools_dedup_primary_only,
        )
        UMI_DEDUP_GENOME = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME
        ch_genomic_dedup_log = UMI_DEDUP_GENOME.out.deduplog
        ch_tsv_edit_distance = UMI_DEDUP_GENOME.out.tsv_edit_distance
        ch_tsv_per_umi = UMI_DEDUP_GENOME.out.tsv_per_umi
        ch_tsv_umi_per_position = UMI_DEDUP_GENOME.out.tsv_umi_per_position
        ch_genome_dedup = UMI_DEDUP_GENOME.out.results.map { r -> [r.id, r.bam, r.bai, r.dedup_log, r.samtools, r.tsv] }
    }

    // Co-ordinate sort, index and run stats on transcriptome BAM. This takes
    // some preparation- we have to coordinate sort the BAM, run the
    // deduplication, then restore name sorting and run a script from umitools
    // to prepare for rsem or salmon

    // 1. Coordinate sort

    BAM_SORT_STATS_SAMTOOLS(
        ch_transcriptome_bam,
        ch_transcript_fasta_fai,
    )
    ch_sorted_transcriptome_bam = BAM_SORT_STATS_SAMTOOLS.out.bam.join(BAM_SORT_STATS_SAMTOOLS.out.index)

    // 2. Transcriptome BAM deduplication
    if (umi_dedup_tool == "umicollapse") {
        BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE_TRANSCRIPTOME(
            ch_sorted_transcriptome_bam
        )
        UMI_DEDUP_TRANSCRIPTOME = BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE_TRANSCRIPTOME
        ch_transcriptomic_dedup_log = UMI_DEDUP_TRANSCRIPTOME.out.dedup_stats
        ch_transcriptome_dedup = UMI_DEDUP_TRANSCRIPTOME.out.results.map { r -> [r.id, r.bam, r.bai, r.dedup_stats, r.samtools, null] }
    }
    else if (umi_dedup_tool == "umitools") {
        BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME(
            ch_sorted_transcriptome_bam,
            umitools_dedup_stats,
            umitools_dedup_primary_only,
        )
        UMI_DEDUP_TRANSCRIPTOME = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME
        ch_transcriptomic_dedup_log = UMI_DEDUP_TRANSCRIPTOME.out.deduplog
        ch_tsv_edit_distance = ch_tsv_edit_distance.mix(UMI_DEDUP_TRANSCRIPTOME.out.tsv_edit_distance)
        ch_tsv_per_umi = ch_tsv_per_umi.mix(UMI_DEDUP_TRANSCRIPTOME.out.tsv_per_umi)
        ch_tsv_umi_per_position = ch_tsv_umi_per_position.mix(UMI_DEDUP_TRANSCRIPTOME.out.tsv_umi_per_position)
        ch_transcriptome_dedup = UMI_DEDUP_TRANSCRIPTOME.out.results.map { r -> [r.id, r.bam, r.bai, r.dedup_log, r.samtools, r.tsv] }
    }

    // 3. Restore name sorting
    SAMTOOLS_SORT(
        UMI_DEDUP_TRANSCRIPTOME.out.bam,
        ch_fasta_fai,
        '',
    )

    // 4. Run prepare_for_rsem.py on paired-end BAM files
    // This fixes paired-end reads in name sorted BAM files
    // See: https://github.com/nf-core/rnaseq/issues/828
    ended_transcriptome_dedup_bam = SAMTOOLS_SORT.out.bam.branch { meta, bam ->
        single_end: meta.single_end
        return [meta, bam]
        paired_end: !meta.single_end
        return [meta, bam]
    }

    UMITOOLS_PREPAREFORRSEM(
        ended_transcriptome_dedup_bam.paired_end.map { meta, bam -> [meta, bam, []] }
    )

    ch_dedup_transcriptome_bam = ended_transcriptome_dedup_bam.single_end.mix(UMITOOLS_PREPAREFORRSEM.out.bam)

    // The transcriptome side is packed into a single element so a remainder
    // join against the genome side yields one null when it is absent (e.g.
    // no transcriptome BAM for HISAT2). Only paired-end samples pass through
    // UMITOOLS_PREPAREFORRSEM.
    ch_transcriptome_results = ch_transcriptome_dedup
        .join(SAMTOOLS_SORT.out.bam.map { meta, bam -> [meta.id, bam] }, by: [0])
        .join(ch_dedup_transcriptome_bam.map { meta, bam -> [meta.id, bam] }, by: [0])
        .join(UMITOOLS_PREPAREFORRSEM.out.bam.map { meta, bam -> [meta.id, bam] }, by: [0], remainder: true)
        .join(UMITOOLS_PREPAREFORRSEM.out.log.map { meta, log -> [meta.id, log] }, by: [0], remainder: true)
        .map { id, dedup_bam, bai, dedup_log, samtools, tsv, sorted_bam, bam, filtered_bam, rsem_log ->
            [
                id,
                record(
                    dedup_log:     dedup_log,
                    rsem_log:      rsem_log,
                    transcriptome: record(
                        bam:              bam,
                        dedup_bam:        dedup_bam,
                        sorted_bam:       sorted_bam,
                        sorted_bam_index: bai,
                        filtered_bam:     filtered_bam,
                        stats:            samtools.stats,
                        flagstat:         samtools.flagstat,
                        idxstats:         samtools.idxstats,
                        tsv:              tsv ? record(edit_distance: tsv.edit_distance, per_umi: tsv.per_umi, umi_per_position: tsv.umi_per_position) : null
                    )
                )
            ]
        }

    ch_results = UMI_DEDUP_GENOME.out.bam
        .map { meta, _bam -> [meta.id, meta] }
        .join(ch_genome_dedup, by: [0])
        .join(ch_transcriptome_results, by: [0], remainder: true)
        .map { id, meta, bam, bai, dedup_log, samtools, tsv, transcriptome ->
            record(
                id:                       id,
                meta:                     meta,
                bam:                      bam,
                bai:                      bai,
                genomic_dedup_log:        dedup_log,
                transcriptomic_dedup_log: transcriptome?.dedup_log,
                prepare_for_rsem_log:     transcriptome?.rsem_log,
                genome:                   record(stats: samtools.stats, flagstat: samtools.flagstat, idxstats: samtools.idxstats),
                transcriptome:            transcriptome?.transcriptome,
                tsv:                      tsv ? record(edit_distance: tsv.edit_distance, per_umi: tsv.per_umi, umi_per_position: tsv.umi_per_position) : null
            ) as UmiDedupBam
        }

    // Collect files useful for MultiQC into one helpful emission. Don't
    // automatically add transcriptome stats- difficult to separate in multiqc
    // without a bit more work

    ch_multiqc_files = ch_genomic_dedup_log
        .mix(UMI_DEDUP_GENOME.out.stats)
        .mix(UMI_DEDUP_GENOME.out.flagstat)
        .mix(UMI_DEDUP_GENOME.out.idxstats)
        .transpose()

    // Genome-side only; transcriptome stats excluded (MultiQC can't
    // disambiguate them from genome stats without extra work).
    ch_per_sample_mqc_bundle = ch_genomic_dedup_log
        .join(UMI_DEDUP_GENOME.out.stats,    remainder: true)
        .join(UMI_DEDUP_GENOME.out.flagstat, remainder: true)
        .join(UMI_DEDUP_GENOME.out.idxstats, remainder: true)
        .map { row -> [row[0], row.drop(1).findAll { f -> f != null }.collectMany { e -> (e instanceof List) ? e : [e] }] }

    emit:
    bam                            = UMI_DEDUP_GENOME.out.bam // channel: [ val(meta), path(bam) ]
    index                          = UMI_DEDUP_GENOME.out.index // channel: [ val(meta), path(bai) ]
    genomic_dedup_log              = ch_genomic_dedup_log // channel: [ val(meta), path(log) ]
    transcriptomic_dedup_log       = ch_transcriptomic_dedup_log // channel: [ val(meta), path(log) ]
    prepare_for_rsem_log           = UMITOOLS_PREPAREFORRSEM.out.log // channel: [ val(meta), path(log) ]
    stats                          = UMI_DEDUP_GENOME.out.stats.mix(UMI_DEDUP_TRANSCRIPTOME.out.stats) // channel: [ val(meta), path(stats)]
    flagstat                       = UMI_DEDUP_GENOME.out.flagstat.mix(UMI_DEDUP_TRANSCRIPTOME.out.flagstat) // channel: [ val(meta), path(flagstat)]
    idxstats                       = UMI_DEDUP_GENOME.out.idxstats.mix(UMI_DEDUP_TRANSCRIPTOME.out.idxstats) // channel: [ val(meta), path(idxstats)]
    genome_stats                   = UMI_DEDUP_GENOME.out.stats // channel: [ val(meta), path(stats) ]
    genome_flagstat                = UMI_DEDUP_GENOME.out.flagstat // channel: [ val(meta), path(flagstat) ]
    genome_idxstats                = UMI_DEDUP_GENOME.out.idxstats // channel: [ val(meta), path(idxstats) ]
    tsv_edit_distance              = ch_tsv_edit_distance // channel: [ val(meta), path(tsv) ]
    tsv_per_umi                    = ch_tsv_per_umi // channel: [ val(meta), path(tsv) ]
    tsv_umi_per_position           = ch_tsv_umi_per_position // channel: [ val(meta), path(tsv) ]
    multiqc_files                  = ch_multiqc_files // channel: [ val(meta), path(file) ]
    transcriptome_bam              = ch_dedup_transcriptome_bam // channel: [ val(meta), path(bam) ] - final output
    transcriptome_dedup_bam        = UMI_DEDUP_TRANSCRIPTOME.out.bam // channel: [ val(meta), path(bam) ] - after dedup, before name sort
    transcriptome_sorted_bam       = SAMTOOLS_SORT.out.bam // channel: [ val(meta), path(bam) ] - name-sorted
    transcriptome_sorted_bam_index = UMI_DEDUP_TRANSCRIPTOME.out.index // channel: [ val(meta), path(index) ] - coordinate-sorted dedup index
    transcriptome_filtered_bam     = UMITOOLS_PREPAREFORRSEM.out.bam // channel: [ val(meta), path(bam) ] - paired-end filtered
    per_sample_mqc_bundle          = ch_per_sample_mqc_bundle // channel: [ val(meta), list(files) ]
    results                        = ch_results // channel: UmiDedupBam
}
