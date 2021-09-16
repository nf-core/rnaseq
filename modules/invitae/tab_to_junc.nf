//
// Convert tab file from STAR to junc format
//

process TAB_TO_JUNC {
    input:
    path tab

    """
    cp $tab > $tab\.junc
    """
}
