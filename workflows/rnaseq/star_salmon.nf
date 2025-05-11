//
// Alignment with STAR and gene/transcript quantification with Salmon
//

include { STAR_ALIGN              } from '../../modules/nf-core/star/align'
include { STAR_ALIGN_IGENOMES     } from '../../modules/local/star_align_igenomes'

include { SAMTOOLS_SORT           } from '../../modules/nf-core/samtools/sort'
include { SAMTOOLS_INDEX          } from '../../modules/nf-core/samtools/index'

include { SAMTOOLS_STATS          } from '../../modules/nf-core/samtools/stats'
include { SAMTOOLS_IDXSTATS       } from '../../modules/nf-core/samtools/idxstats'
include { SAMTOOLS_FLAGSTAT       } from '../../modules/nf-core/samtools/flagstat'

include { UMITOOLS_DEDUP          } from '../../modules/nf-core/umitools/dedup'
include { UMITOOLS_PREPAREFORRSEM as UMITOOLS_PREPAREFORSALMON } from '../../modules/nf-core/umitools/prepareforrsem'

include { SALMON_QUANT            } from '../../modules/nf-core/salmon/quant'
include { CUSTOM_TX2GENE          } from '../../modules/nf-core/custom/tx2gene'
include { TXIMETA_TXIMPORT        } from '../../modules/nf-core/tximeta/tximport'

include {
    SUMMARIZEDEXPERIMENT_SUMMARIZEDEXPERIMENT as SE_GENE ;
    SUMMARIZEDEXPERIMENT_SUMMARIZEDEXPERIMENT as SE_GENE_LENGTH_SCALED ;
    SUMMARIZEDEXPERIMENT_SUMMARIZEDEXPERIMENT as SE_GENE_SCALED ;
    SUMMARIZEDEXPERIMENT_SUMMARIZEDEXPERIMENT as SE_TRANSCRIPT
} from '../../modules/nf-core/summarizedexperiment/summarizedexperiment'


workflow STAR_SALMON {
    take:
    reads               // channel: [ val(meta), [ reads ] ]
    index               // channel: [ val(meta), [ index ] ]
    gtf                 // channel: [ val(meta), [ gtf ] ]
    star_ignore_sjdbgtf // boolean: when using pre-built STAR indices do not re-extract and use splice junctions from the GTF file
    seq_platform        // string : sequencing platform
    seq_center          // string : sequencing center
    is_aws_igenome      // boolean: whether the genome files are from AWS iGenomes
    fasta               // channel: /path/to/fasta
    transcript_fasta          // channel: /path/to/transcript.fasta
    gtf_id_attribute          //     val: GTF gene ID attribute
    gtf_extra_attribute       //     val: GTF alternative gene attribute (e.g. gene_name)
    lib_type                  //     val: String to override Salmon library type

    main:
    bam = reads.map { meta, fastq ->
        def star = is_aws_igenome
            ? STAR_ALIGN_IGENOMES ( meta, fastq, index, gtf, star_ignore_sjdbgtf, seq_platform, seq_center )
            : STAR_ALIGN ( meta, fastq, index, gtf, star_ignore_sjdbgtf, seq_platform, seq_center )

        def bam_sort = SAMTOOLS_SORT ( meta, star.bam, fasta )
        def bam_index = SAMTOOLS_INDEX ( meta, bam_sort.bam )
        def bam_bai_csi = bam_index.bai ?: bam_index.csi
        def bam_stats = [
            SAMTOOLS_STATS ( meta, bam_index.bam, bam_bai_csi, fasta ),
            SAMTOOLS_FLAGSTAT ( meta, bam_index.bam, bam_bai_csi ),
            SAMTOOLS_IDXSTATS ( meta, bam_index.bam, bam_bai_csi ),
        ]

        def bam_tx = star.bam_transcript

        //
        // SUBWORKFLOW: Remove duplicate reads from BAM file based on UMIs
        //
        if (params.with_umi) {
            // Deduplicate genome BAM file before downstream analysis
            def umi_dedup = UMITOOLS_DEDUP ( meta, bam_index.bam, bam_index.bai, params.umitools_dedup_stats )
            def umi_index = SAMTOOLS_INDEX ( meta, umi_dedup.bam )
            def umi_bai_csi = umi_index.bai ?: umi_index.csi
            def umi_stats = [
                SAMTOOLS_STATS ( meta, umi_index.bam, umi_bai_csi, fasta ),
                SAMTOOLS_FLAGSTAT ( meta, umi_index.bam, umi_bai_csi ),
                SAMTOOLS_IDXSTATS ( meta, umi_index.bam, umi_bai_csi ),
            ]

            // TODO: add umi results to record type
            [ meta, umi_dedup, umi_index, umi_stats ]

            // Co-ordinate sort, index and run stats on transcriptome BAM
            def bam_tx_sort = SAMTOOLS_SORT ( meta, bam_tx, transcript_fasta )
            def bam_tx_index = SAMTOOLS_INDEX ( meta, bam_tx_sort.bam )

            // Deduplicate transcriptome BAM file before read counting with Salmon
            def umi_tx_dedup = UMITOOLS_DEDUP ( meta, bam_tx_index.bam, bam_tx_index.bai, params.umitools_dedup_stats )

            // Name sort BAM before passing to Salmon
            // Fix paired-end reads in name sorted BAM file
            // See: https://github.com/nf-core/rnaseq/issues/828
            // Only run prepare_for_rsem.py on paired-end BAM files
            def umi_tx_sort = SAMTOOLS_SORT ( meta, umi_tx_dedup.bam, fasta )

            bam_tx = meta.single_end
                ? umi_tx_sort.bam
                : UMITOOLS_PREPAREFORSALMON ( meta, umi_tx_sort.bam )
        }

        //
        // SUBWORKFLOW: Count reads from BAM alignments using Salmon
        //
        def quant = SALMON_QUANT ( meta, fastq, null, gtf, transcript_fasta, true, lib_type )

        // TODO: StarSalmonSample record type
        [ meta, star, bam_sort, bam_index, bam_stats, quant ]
    }

    tx2gene = CUSTOM_TX2GENE (
        [:],
        gtf,
        bam.map { sample -> sample.quant }.collect(),
        'salmon',
        gtf_id_attribute,
        gtf_extra_attribute
    )

    tximport = TXIMETA_TXIMPORT (
        ['id': 'all_samples'],
        bam.map { sample -> sample.quant }.collect(),
        tx2gene,
        'salmon'
    )

    se_gene = SE_GENE (
        ['id': 'all_samples'],
        [ tximport.counts_gene, tximport.tpm_gene ],
        tx2gene,
        null // samplesheet
    )

    se_gene_length_scaled = SE_GENE_LENGTH_SCALED (
        ['id': 'all_samples'],
        [ tximport.counts_gene_length_scaled, tximport.tpm_gene ],
        tx2gene,
        null // samplesheet
    )

    se_gene_scaled = SE_GENE_SCALED (
        ['id': 'all_samples'],
        [ tximport.counts_gene_scaled, tximport.tpm_gene ],
        tx2gene,
        null // samplesheet
    )

    se_transcript = SE_TRANSCRIPT (
        ['id': 'all_samples'],
        [ tximport.counts_transcript, tximport.tpm_transcript ],
        tx2gene,
        null // samplesheet
    )

    // TODO: StarSalmonCounts record type
    counts = [
        tpm_gene:                      tximport.tpm_gene,
        counts_gene:                   tximport.counts_gene,
        lengths_gene:                  tximport.lengths_gene,
        counts_gene_length_scaled:     tximport.counts_gene_length_scaled,
        counts_gene_scaled:            tximport.counts_gene_scaled,
        tpm_transcript:                tximport.tpm_transcript,
        counts_transcript:             tximport.counts_transcript,
        lengths_transcript:            tximport.lengths_transcript,
        merged_gene_rds:               se_gene.rds,
        merged_gene_rds_length_scaled: se_gene_length_scaled.rds,
        merged_gene_rds_scaled:        se_gene_scaled.rds,
        merged_transcript_rds:         se_transcript.rds,
    ]

    emit:
    bam
    counts

}
