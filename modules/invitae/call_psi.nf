//
// 1. Convert tab file from STAR to junc format
// 2. Invoke PSI caller on this junc file, with respect to reference PON
//

process CALL_PSI {
    input:
    tuple val(meta), path(tab)

    output:
    tuple val(meta), path('*.junc'), emit: junc
    tuple val(meta), path('*.bed'), emit: bed

    """
    date > date.txt
    awk -F '\t' -v OFS='\t' '{if(\$4==0) strand="."; else if(\$4==1) strand="+"; else strand="-"; print \$1, \$2, \$3, ".", \$7, strand;}' $tab > ${tab}.junc
    date >> date.txt
    source /locus/data/dev_analysis/pipe-dev-2/opt/locus-setup/v2/setup.sh nkampshughes/git/pipe_splice_research/locus-pipe
    date >> date.txt
    /locus/home/nkampshughes/git/pipe_splice_research/locus-pipe/bin/germline-rna/psi-junction-caller \\
        --annotated-junctions /locus/home/rmartin/UCSC/transcriptome_model/UCSC_ANNOTATEDJUNCTIONS_INCREMENT_EXON_START.bed \\
        --input-pon /locus/home/rmartin/UCSC/transcriptome_model/UCSC_PON_for_SQ9027_INCREMENT_EXON_START_chrALL.bed \\
        --targets /locus/home/rmartin/UCSC/transcriptome_model/UCSC_TARGETS_INCREMENT_EXON_START.bed \\
        --output-junctions ${tab}.bed \\
        --input-junctions ${tab}.junc
    date >> date.txt
    """
}
