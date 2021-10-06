//
// 1. Convert tab file from STAR to junc format
// 2. Invoke PSI caller on this junc file, with respect to reference PON
//

process CALL_PSI {
    input:
    tuple val(meta), path(tab)
    path template_script

//legacy
//    output:
//    tuple val(meta), path('*.junc'), emit: junc
//    tuple val(meta), path('*.bed'), emit: bed

    """
    date > date.txt
    awk -F '\t' -v OFS='\t' '{if(\$4==0) strand="."; else if(\$4==1) strand="+"; else strand="-"; print \$1, \$2, \$3, ".", \$7, strand;}' $tab > ${tab}.junc
    date >> date.txt
    awk -v junc='${tab}.junc' \\
        '{gsub("JUNC_REPLACE_ME",junc); print \$0;}' \\
        $template_script > ${template_script}_ready.sh
    date >> date.txt
    qsub ${template_script}_ready.sh
    """
}
