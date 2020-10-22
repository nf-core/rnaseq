// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

params.options = [:]

/*
 * Fetch SRA / ENA / GEO run information via the ENA API
 */
process SRA_IDS_TO_RUNINFO {
    tag "$id"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "conda-forge::requests=2.24.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/requests:2.24.0"
    } else {
        // container "quay.io/biocontainers/requests:2.24.0"
        // container "quay.io/biocontainers/requests:2.19.1"
        // Python and requests
        // container "quay.io/biocontainers/mulled-v2-ffdffc678ef7e057a54c6e2a990ebda211c39d9c:62ed692956f96d8218414426148c868f259f179c-0"
        container "nfcore/viralrecon:1.1.0"
    }
    
    input:
    val id
    
    output:
    path "*.tsv", emit: tsv
    
    script:
    """
    echo $id > id.txt
    sra_ids_to_runinfo.py id.txt ${id}_run_info.tsv
    """
}
