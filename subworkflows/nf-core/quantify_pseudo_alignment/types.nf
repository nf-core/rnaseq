// Documentation only: nothing casts a record(...) to these types (nextflow-io/nextflow#7680 corrupts remote Path fields on cast).
include { QuantMerged } from '../quant_tximport_summarizedexperiment/types'

record PseudoQuantSample {
    id:           String
    meta:         Map
    quant_dir:    Path
    json_info:    Path?
    log:          Path?
    quant_merged: QuantMerged
}
