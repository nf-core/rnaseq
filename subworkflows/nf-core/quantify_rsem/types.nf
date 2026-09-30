// Documentation only: nothing casts a record(...) to these types (nextflow-io/nextflow#7680 corrupts remote Path fields on cast).
include { QuantMerged } from '../quant_tximport_summarizedexperiment/types'

record Reads {
    id:    String
    meta:  Map
    reads: List<Path>
}

record RsemMerge {
    counts_gene:       Path
    tpm_gene:          Path
    counts_transcript: Path
    tpm_transcript:    Path
    genes_long:        Path
    isoforms_long:     Path
}

record RsemQuantSample {
    id:                String
    meta:              Map
    counts_gene:       Path
    counts_transcript: Path
    stat:              Path
    log:               Path?
    bam_star:          Path?
    bam_genome:        Path?
    bam_transcript:    Path?
}

// One row per sample (sample set, merged and rsem_merge null) plus the merged rows (sample null):
// one per sample under skip_merge, a single 'all_samples' row otherwise. rsem_merge is null under skip_merge.
record RsemQuantified {
    id:         String
    meta:       Map
    sample:     RsemQuantSample?
    merged:     QuantMerged?
    rsem_merge: RsemMerge?
}
