// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

process FEATURECOUNTS_MERGE_COUNTS {
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container "biocontainers/biocontainers:v1.2.0_cv1"

    conda (params.conda ? "conda-forge::sed=4.7" : null)
    
    input:
    path ('counts/*')
    val  options

    output:
    path "*.counts.tsv", emit: counts

    script:
    """
    mkdir tmp_counts
    cut -f 1,7 `ls ./counts/* | head -n 1` | grep -v "^#" | tail -n+1 > ids.tsv
    for fileid in `ls ./counts/*`; do
        samplename=`basename \$fileid | sed s/\\.featureCounts.txt\$//g`
        echo \$samplename > tmp_counts/\$samplename.featureCounts.txt
        grep -v "^#" \${fileid} | cut -f 8 | tail -n+2 >> tmp_counts/\$samplename.featureCounts.txt
    done

    paste ids.tsv tmp_counts/* > featurecounts.merged.counts.tsv
    """
}
