//
// Convert tab file from STAR to junc format
//

process JUNC_FROM_TAB {
    input:
    tuple val(meta), path(tab)

    output:
    tuple val(meta), path('*.junc'), emit: junc

    """
    date > date.txt
    awk -F '\t' -v OFS='\t' '{if(\$4==0) strand="."; else if(\$4==1) strand="+"; else strand="-"; print \$1, \$2, \$3, ".", \$7, strand;}' $tab > ${tab}.junc
    """
}
