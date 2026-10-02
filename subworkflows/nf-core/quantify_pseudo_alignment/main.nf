nextflow.enable.types = true

//
// Pseudoalignment and quantification with Salmon or Kallisto
//

include { SALMON_QUANT     } from '../../../modules/nf-core/salmon/quant'
include { KALLISTO_QUANT   } from '../../../modules/nf-core/kallisto/quant'

include { QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT } from '../quant_tximport_summarizedexperiment'
include { ReadsInput; QuantMerged } from '../../../modules/nf-core/types'
include { KallistoQuantSample } from '../../../modules/nf-core/kallisto/quant/main'
include { SalmonQuantSample } from '../../../modules/nf-core/salmon/quant/main'

//
// Salmon quant options: the library type is the given one if it is valid, else derived from the sample strandedness
//
def salmonQuantArgs(libtype: String?, extra_args: String?, meta: Map) -> String {
    def valid_libtypes = [ 'A', 'U', 'SF', 'SR', 'IS', 'IU', 'ISF', 'ISR', 'OS', 'OU', 'OSF', 'OSR', 'MS', 'MU', 'MSF', 'MSR' ]
    def lib_type = 'A'
    if (libtype) {
        if (valid_libtypes.contains(libtype)) {
            lib_type = libtype
        }
        else {
            log.info("[Salmon Quant] Invalid library type specified '--libType=${libtype}', defaulting to auto-detection with '--libType=A'.")
        }
    }
    else {
        def strandedness = meta.strandedness
        lib_type = meta.single_end
            ? (strandedness == 'forward' ? 'SF' : strandedness == 'reverse' ? 'SR' : 'U')
            : (strandedness == 'forward' ? 'ISF' : strandedness == 'reverse' ? 'ISR' : 'IU')
    }
    return [ "--libType=${lib_type}", extra_args ?: '' ].join(' ').trim()
}

//
// Kallisto quant options: the strandedness of the sample unless the extra options set one
//
def kallistoQuantArgs(extra_args: String?, meta: Map) -> String {
    def extra = extra_args ?: ''
    def strandedness = (extra.contains('--fr-stranded') || extra.contains('--rf-stranded')) ? '' :
        meta.strandedness == 'forward' ? '--fr-stranded' :
        meta.strandedness == 'reverse' ? '--rf-stranded' : ''
    return [ strandedness, extra ].join(' ').trim()
}

// The Salmon and Kallisto quant options, and the prefix of the merged SummarizedExperiment files
record QuantPseudoArgs {
    salmon_quant:   String?
    kallisto_quant: String?
    se_prefix:      String?
}

workflow QUANTIFY_PSEUDO_ALIGNMENT {
    take:
    samplesheet: Value<Path>
    ch_samples: Channel<ReadsInput>
    index: Value<Path?>
    transcript_fasta: Value<Path?>
    gtf: Value<Path>
    gtf_id_attribute: String // GTF gene ID attribute
    gtf_extra_attribute: String // GTF alternative gene attribute (e.g. gene_name)
    pseudo_aligner: String // kallisto or salmon
    kallisto_quant_fraglen: Integer? // Estimated fragment length required by Kallisto in single-end mode
    kallisto_quant_fraglen_sd: Integer? // Estimated standard error for fragment length required by Kallisto in single-end mode
    skip_merge: Boolean // skip cross-sample merging, run tximport per-sample
    salmon_libtype: String? // Salmon library type; derived from the sample strandedness when null
    tool_args: QuantPseudoArgs // extra Salmon and Kallisto quant options (salmon_quant, kallisto_quant) and the prefix of the merged SummarizedExperiment files (se_prefix, default: the sample id)

    main:

    //
    // Quantify and merge counts across samples
    //
    def ch_salmon: Channel<SalmonQuantSample>     = channel.empty()
    def ch_kallisto: Channel<KallistoQuantSample> = channel.empty()
    if (pseudo_aligner == 'salmon') {
        ch_salmon = SALMON_QUANT(
            ch_samples.map { r -> r + record(args: salmonQuantArgs(salmon_libtype, tool_args.salmon_quant, r.meta)) },
            index,
            gtf,
            transcript_fasta
        )
    } else {
        ch_kallisto = KALLISTO_QUANT(
            ch_samples.map { r -> r + record(args: kallistoQuantArgs(tool_args.kallisto_quant, r.meta)) },
            index,
            gtf,
            null,
            kallisto_quant_fraglen,
            kallisto_quant_fraglen_sd
        )
    }
    def ch_quant_dirs = ch_salmon.mix(ch_kallisto).map { r -> record(id: r.id, meta: r.meta, quants: [ r.quant_dir ]) }

    //
    // Post-process quantifications with tximport and SummarizedExperiment
    //
    def ch_quant_merged: Channel<QuantMerged> = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT(
        samplesheet,
        ch_quant_dirs,
        gtf,
        gtf_id_attribute,
        gtf_extra_attribute,
        pseudo_aligner,
        skip_merge,
        tool_args
    )

    emit:
    salmon:   Channel<SalmonQuantSample>   = ch_salmon        // per sample, when pseudo_aligner is salmon
    kallisto: Channel<KallistoQuantSample> = ch_kallisto      // per sample, when pseudo_aligner is kallisto
    merged:   Channel<QuantMerged>         = ch_quant_merged  // one row per sample under skip_merge, a single 'all_samples' row otherwise
}
