include { QuantMerged } from '../quant_tximport_summarizedexperiment/types'

record PseudoQuantSample {
    id:           String
    meta:         Map
    quant_dir:    Path
    json_info:    Path?
    log:          Path?
    quant_merged: QuantMerged
}
