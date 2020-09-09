// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process SUBREAD_FEATURECOUNTS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "quay.io/biocontainers/subread:2.0.1--hed695b0_0"
    //container "https://depot.galaxyproject.org/singularity/subread:2.0.1--hed695b0_0"

    conda (params.conda ? "bioconda::subread=2.0.1" : null)

    input:
    tuple val(meta), path(bams), path(annotation)
    val options

    output:
    tuple val(meta), path("*featureCounts.txt"), emit: txt
    tuple val(meta), path("*featureCounts.txt.summary"), emit: summary
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    def prefix   = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"

    def strandedness = 0
    if (meta.strandedness == 'forward') {
        strandedness = 1
    } else if (meta.strandedness == 'reverse') {
        strandedness = 2
    }
    def paired_end = meta.single_end ? '' : '-p'
    """
    featureCounts \\
        $ioptions.args \\
        $paired_end \\
        -T $task.cpus \\
        -a $annotation \\
        -s $strandedness \\
        -o ${prefix}.featureCounts.txt \\
        ${bams.join(' ')}

    echo \$(featureCounts -v 2>&1) | sed -e "s/featureCounts v//g" > ${software}.version.txt
    """
}

// input:
// path bam from bam_featurecounts
// path gtf from ch_gtf
// path biotypes_header from ch_biotypes_header

// output:
// path "${bam.baseName}_gene.featureCounts.txt" into geneCounts,
//                                                    featureCounts_to_merge
// path "${bam.baseName}_gene.featureCounts.txt.summary" into featureCounts_logs
// path "${bam.baseName}_biotype_counts*mqc.{txt,tsv}" optional true into featureCounts_biotype

// script:
// // Try to get real sample name
// sample_name = bam.baseName - 'Aligned.sortedByCoord.out' - '_subsamp.sorted'
// biotype_qc = params.skip_biotype_qc ? '' : "featureCounts -a $gtf -g $biotype -t ${params.fc_count_type} -o ${bam.baseName}_biotype.featureCounts.txt -p -s $featureCounts_direction $bam"
// mod_biotype = params.skip_biotype_qc ? '' : "cut -f 1,7 ${bam.baseName}_biotype.featureCounts.txt | tail -n +3 | cat $biotypes_header - >> ${bam.baseName}_biotype_counts_mqc.txt && mqc_features_stat.py ${bam.baseName}_biotype_counts_mqc.txt -s $sample_name -f rRNA -o ${bam.baseName}_biotype_counts_gs_mqc.tsv"
// """
// featureCounts \\
//     -a $gtf \\
//     -g $params.fc_group_features \\
//     -t $params.fc_count_type \\
//     -o ${bam.baseName}_gene.featureCounts.txt \\
//     $extraAttributes \\
//     -p \\
//     -s $featureCounts_direction \\
//     $bam
// $biotype_qc
// $mod_biotype
// """
