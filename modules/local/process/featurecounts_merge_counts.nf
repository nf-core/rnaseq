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
    path "*.txt", emit: counts

    script:
    """
    mkdir tmp_counts
    cut -f 1,7 ${counts[0]} | grep -v "^#" | tail -n+1 > ids.tsv
    for fileid in $counts; do
        basename \$fileid | sed s/\\.featureCounts.txt\$//g > tmp_counts/\$fileid
        grep -v "^#" \${fileid} | cut -f 8 | tail -n+2 >> tmp_counts/\$fileid
    done

    paste ids.tsv tmp_counts/* > featurecounts.merged.counts.txt
    """
}
