include { BCFTOOLS_MPILEUP     } from '../../../modules/nf-core/bcftools/mpileup/main'
include { NGSCHECKMATE_NCM } from '../../../modules/nf-core/ngscheckmate/ncm/main.nf'

workflow BAM_NGSCHECKMATE {

    take:
    ch_bam              // channel: [ val(meta), bam ]
    ch_bed              // channel: [ bed ]
    ch_fasta            // channel: [ fasta ]

    main:

    ch_versions = Channel.empty()

    ch_bam_bed  = ch_bam.combine(ch_bed)

    BCFTOOLS_MPILEUP (ch_bam_bed, ch_fasta.collect(), false)
    ch_versions = ch_versions.mix(BCFTOOLS_MPILEUP.out.versions)

    BCFTOOLS_MPILEUP
    .out
    .vcf
    .map{meta, vcf -> vcf}
    .collect()
    .set {ch_flat_vcfs}

    NGSCHECKMATE_NCM (ch_flat_vcfs, ch_bed, ch_fasta)
    ch_versions = ch_versions.mix(NGSCHECKMATE_NCM.out.versions)

    emit:
    pdf          = NGSCHECKMATE_NCM.out.pdf          // channel: [ pdf ]
    corr_matrix  = NGSCHECKMATE_NCM.out.corr_matrix  // channel: [ corr_matrix ]
    matched      = NGSCHECKMATE_NCM.out.matched      // channel: [ matched ]
    all          = NGSCHECKMATE_NCM.out.all          // channel: [ all ]
    vcf          = BCFTOOLS_MPILEUP.out.vcf          // channel: [ meta, vcf ]
    versions     = ch_versions                       // channel: [ versions.yml ]

}

