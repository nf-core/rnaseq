include { BCFTOOLS_MPILEUP } from '../../../modules/nf-core/bcftools/mpileup/main'
include { NGSCHECKMATE_NCM } from '../../../modules/nf-core/ngscheckmate/ncm/main'

workflow BAM_NGSCHECKMATE {

    take:
    ch_input            // channel: [ val(meta1), bam/cram ]
    ch_snp_bed          // channel: [ val(meta2), bed ]
    ch_fasta            // channel: [ val(meta3), fasta ]

    main:

    ch_versions = Channel.empty()

    ch_input_bed = ch_input.combine(ch_snp_bed.collect())
                        // do something to combine the metas?
                        .map{ input_meta, input_file, bed_meta, bed_file ->
                            [input_meta, input_file, bed_file]
                        }

    BCFTOOLS_MPILEUP (ch_input_bed, ch_fasta.collect(), false)
    ch_versions = ch_versions.mix(BCFTOOLS_MPILEUP.out.versions)

    BCFTOOLS_MPILEUP
    .out
    .vcf
    .map{meta, vcf -> vcf}    // discard individual metas
    .collect()                // group into one channel
    .map{files -> [files]}    // make the channel into [vcf1, vcf2, ...]
    .set {ch_collected_vcfs}

    ch_snp_bed
    .map{meta, bed -> meta} // use the snp_bed file meta as the meta for the merged channel
    .combine(ch_collected_vcfs) // add the vcf files after the meta, now looks like [meta, [vcf1, vcf2, ... ] ]
    .set {ch_vcfs}

    NGSCHECKMATE_NCM (ch_vcfs, ch_snp_bed, ch_fasta)
    ch_versions = ch_versions.mix(NGSCHECKMATE_NCM.out.versions)

    emit:
    corr_matrix  = NGSCHECKMATE_NCM.out.corr_matrix  // channel: [ meta, corr_matrix ]
    matched      = NGSCHECKMATE_NCM.out.matched      // channel: [ meta, matched ]
    all          = NGSCHECKMATE_NCM.out.all          // channel: [ meta, all ]
    vcf          = BCFTOOLS_MPILEUP.out.vcf          // channel: [ meta, vcf ]
    pdf          = NGSCHECKMATE_NCM.out.pdf          // channel: [ meta, pdf ]
    versions     = ch_versions                       // channel: [ versions.yml ]

}

