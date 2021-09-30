//
// 1. Convert tab file from STAR to junc format
// 2. Invoke PSI caller on this junc file, with respect to reference PON
//

include { CALL_PSI        } from '../../modules/invitae/call_psi'

workflow PSI_CALLIN {
    take:
    ch_tab // channel: [ val(meta), [ tab ] ]
    ch_script // file: /path/to/source.sh
    ch_caller // file: /path/to/caller executable
    ch_annotations // file: /path/to/annotated-junctions.bed
    ch_pon // file: /path/to/pon.bed
    ch_targets // file: /path/to/targets.bed

    main:
    CALL_PSI ( ch_tab, ch_script, ch_caller, ch_annotations, ch_pon, ch_targets )

    emit:
    junc         = CALL_PSI.out.junc
    bed         = CALL_PSI.out.bed
}
