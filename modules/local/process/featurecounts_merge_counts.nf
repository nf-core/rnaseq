// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

params.options = [:]

process FEATURECOUNTS_MERGE_COUNTS {
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    if (workflow.containerEngine == 'singularity' || !params.pull_docker_container) {
        container "biocontainers/biocontainers:v1.2.0_cv1"
    } else {
        container "biocontainers/biocontainers:v1.2.0_cv1"
    }
    
    input:
    path ('counts/*')
    
    output:
    path "*.counts.tsv", emit: counts

    script:
    """
    mkdir -p tmp/counts
    cut -f 1,7 `ls ./counts/* | head -n 1` | grep -v "^#" | tail -n+1 > ids.tsv
    for fileid in `ls ./counts/*`; do
        samplename=`basename \$fileid | sed s/\\.featureCounts.txt\$//g`
        echo \$samplename > tmp/counts/\$samplename.featureCounts.txt
        grep -v "^#" \${fileid} | cut -f 8 | tail -n+2 >> tmp/counts/\$samplename.featureCounts.txt
    done

    paste ids.tsv tmp/counts/* > featurecounts.merged.counts.tsv
    """
}
