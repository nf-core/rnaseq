//
// Convert tab file from STAR to junc format
//

include { JUNC_FROM_TAB        } from '../../modules/invitae/junc_from_tab'

workflow TAB_TO_JUNC {
    take:
    tab

    main:
    JUNC_FROM_TAB ( tab )

    emit:
    junc         = JUNC_FROM_TAB.out.junc
}
