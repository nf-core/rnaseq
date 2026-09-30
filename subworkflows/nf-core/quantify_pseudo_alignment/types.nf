// Documentation only: nothing casts a record(...) to these types (nextflow-io/nextflow#7680 corrupts remote Path fields on cast).
include { QuantMerged } from '../quant_tximport_summarizedexperiment/types'

record Reads {
    id:    String
    meta:  Map
    reads: List<Path>
}

// multiqc is the file MultiQC parses: the quant directory for Salmon, the log for Kallisto.
record PseudoQuantSample {
    id:        String
    meta:      Map
    quant_dir: Path
    json_info: Path?
    log:       Path?
    multiqc:   Path
}

// One row per sample (sample set, merged null) plus the merged rows (sample null):
// one per sample under skip_merge, a single 'all_samples' row otherwise.
record PseudoQuantified {
    id:     String
    meta:   Map
    sample: PseudoQuantSample?
    merged: QuantMerged?
}
