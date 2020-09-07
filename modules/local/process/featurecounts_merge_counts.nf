// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

process FEATURECOUNTS_MERGE_COUNTS {
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.conda ? "${baseDir}/environment.yml" : null)

    input:
    path counts
    val  options

    output:
    path ".txt", emit: tsv

    script:
    // Redirection (the `<()`) for the win!
    // Geneid in 1st column and gene_name in 7th
    ids = "<(tail -n +2 ${counts[0]} | cut -f1,7 )"
    // Remove first line and take third column
    files = counts.collect { filename -> "<(tail -n +2 ${filename} | sed 's:.sorted.bam::' | cut -f8)" }.join(" ")
    """
    paste $ids $files > featurecounts.merged.counts.tsv
    """
}
