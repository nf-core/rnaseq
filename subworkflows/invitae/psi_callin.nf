//
// 1. Convert tab file from STAR to junc format
// 2. Invoke PSI caller on this junc file, with respect to reference PON
//

include { CALL_PSI        } from '../../modules/invitae/call_psi'

workflow PSI_CALLIN {
    take:
    ch_tab // channel: [ val(meta), [ tab ] ]

    main:
    CALL_PSI ( ch_tab )

    emit:
    junc         = CALL_PSI.out.junc
    bed         = CALL_PSI.out.bed
}
