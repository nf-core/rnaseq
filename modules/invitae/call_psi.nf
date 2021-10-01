//
// 1. Convert tab file from STAR to junc format
// 2. Invoke PSI caller on this junc file, with respect to reference PON
//

process CALL_PSI {
    input:
    tuple val(meta), path(tab)
    path venv
    path caller
    path annotations
    path pon
    path targets

    output:
    tuple val(meta), path('*.junc'), emit: junc
    tuple val(meta), path('*.bed'), emit: bed

    // previous attempts
    // source /locus/data/dev_analysis/pipe-dev-2/opt/locus-setup/v2/setup.sh nkampshughes/git/pipe_splice_research/locus-pipe
    // this failed with error code 1, suggesting it cannot find the script ...
    // as an attempt at a solution, i copied the script components in a hardcoded format below,
    // where setup.bash is a root directory symbolic link to that hardcoded destination
    // similarly, replacing this with something that can be accessed symbolically:
    // /locus/home/nkampshughes/git/pipe_splice_research/locus-pipe/bin/germline-rna/psi-junction-caller

    """
    date > date.txt
    awk -F '\t' -v OFS='\t' '{if(\$4==0) strand="."; else if(\$4==1) strand="+"; else strand="-"; print \$1, \$2, \$3, ".", \$7, strand;}' $tab > ${tab}.junc
    date >> date.txt
    source $venv
    date >> date.txt
    python $caller \\
        --annotated-junctions $annotations \\
        --input-pon $pon \\
        --targets $targets \\
        --output-junctions ${tab}.bed \\
        --input-junctions ${tab}.junc
    date >> date.txt
    """
}
