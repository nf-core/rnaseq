include { DUPRADAR              } from '../../../modules/nf-core/dupradar'
include { PRESEQ_LCEXTRAP       } from '../../../modules/nf-core/preseq/lcextrap'
include { QUALIMAP_RNASEQ       } from '../../../modules/nf-core/qualimap/rnaseq'
include { SUBREAD_FEATURECOUNTS } from '../../../modules/nf-core/subread/featurecounts'
include { MULTIQC_CUSTOM_BIOTYPE } from '../../../modules/local/multiqc_custom_biotype'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_QUALIMAP } from '../../../modules/nf-core/samtools/sort'
include { BAM_RSEQC             } from '../../nf-core/bam_rseqc'
include { biotypeInGtf          } from '../utils_nfcore_rnaseq_pipeline'

workflow BAM_POST_ALIGNMENT_QC {

    take:
    ch_genome_bam           // [ meta, bam ]
    ch_genome_bam_index     // [ meta, bai ]
    ch_gtf                  // path(gtf)
    ch_gene_bed             // path(gene.bed)
    ch_fasta_fai            // [ meta, fasta, fai ] for samtools sort
    ch_biotypes_header      // path(biotypes_header.txt)
    skip_preseq             // val(bool)
    skip_biotype_qc         // val(bool)
    skip_qualimap           // val(bool)
    skip_dupradar           // val(bool)
    skip_rseqc              // val(bool)
    biotype                 // val(string) - e.g. "gene_type" or "gene_biotype"
    rseqc_modules           // val(list) - pre-parsed list of RSeQC module names

    main:
    ch_multiqc_files       = channel.empty()
    ch_inferexperiment_txt = channel.empty()
    ch_versions            = channel.empty()

    //
    // MODULE: Run Preseq
    //
    if (!skip_preseq) {
        PRESEQ_LCEXTRAP (
            ch_genome_bam
        )
        ch_multiqc_files = ch_multiqc_files.mix(PRESEQ_LCEXTRAP.out.lc_extrap)
    }

    //
    // MODULE: Feature biotype QC using featureCounts
    //
    if (!skip_biotype_qc && biotype) {

        ch_gtf
            .map { gtf -> biotypeInGtf(gtf, biotype) }
            .set { biotype_in_gtf }

        // Prevent any samples from running if GTF file doesn't have a valid biotype
        ch_genome_bam
            .combine(ch_gtf)
            .combine(biotype_in_gtf)
            .filter { meta, bam, gtf, biotype_ok -> biotype_ok }
            .map { meta, bam, gtf, _biotype_ok -> [ meta, bam, gtf ] }
            .set { ch_featurecounts }

        SUBREAD_FEATURECOUNTS (
            ch_featurecounts
        )

        MULTIQC_CUSTOM_BIOTYPE (
            SUBREAD_FEATURECOUNTS.out.counts,
            ch_biotypes_header
        )
        ch_multiqc_files = ch_multiqc_files.mix(MULTIQC_CUSTOM_BIOTYPE.out.tsv)
    }

    //
    // MODULE: Qualimap
    //
    if (!skip_qualimap) {
        // Sort BAM by name for qualimap (performance optimization)
        SAMTOOLS_SORT_QUALIMAP (
            ch_genome_bam,
            ch_fasta_fai,
            ''
        )

        QUALIMAP_RNASEQ (
            SAMTOOLS_SORT_QUALIMAP.out.bam,
            ch_gtf.map { item -> [ [:], item ] }
        )
        ch_multiqc_files = ch_multiqc_files.mix(QUALIMAP_RNASEQ.out.results)
    }

    //
    // MODULE: dupRadar
    //
    if (!skip_dupradar) {
        DUPRADAR (
            ch_genome_bam,
            ch_gtf.map { item -> [ [:], item ] }
        )
        ch_multiqc_files = ch_multiqc_files.mix(DUPRADAR.out.multiqc)
        ch_versions = ch_versions.mix(DUPRADAR.out.versions.first())
    }

    //
    // SUBWORKFLOW: RSeQC
    //
    if (!skip_rseqc && rseqc_modules.size() > 0) {
        BAM_RSEQC (
            ch_genome_bam.join(ch_genome_bam_index, by: [0]).map { meta, bam, bai -> [ meta, [ bam, bai ] ] },
            ch_gene_bed,
            rseqc_modules
        )
        ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.bamstat_txt)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.inferexperiment_txt)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.innerdistance_freq)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.junctionannotation_log)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.junctionsaturation_rscript)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.readdistribution_txt)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.readduplication_pos_xls)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_RSEQC.out.tin_txt)
        ch_inferexperiment_txt = BAM_RSEQC.out.inferexperiment_txt
    }

    emit:
    multiqc_files       = ch_multiqc_files
    inferexperiment_txt = ch_inferexperiment_txt
    versions            = ch_versions
}
