// Import generic module functions
include { getSoftwareName; getProcessName } from "$projectDir/lib/functions"

def VERSION = '2.2.0'

process HISAT2_EXTRACTSPLICESITES {
    tag "$gtf"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::hisat2=2.2.1' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/hisat2:2.2.1--h1b792b2_3"
    } else {
        container "quay.io/biocontainers/hisat2:2.2.1--h1b792b2_3"
    }

    input:
    path gtf

    output:
    path "*.splice_sites.txt", emit: txt
    path "versions.yml"      , emit: versions

    script:
    """
    hisat2_extract_splice_sites.py $gtf > ${gtf.baseName}.splice_sites.txt
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo $VERSION)
    END_VERSIONS
    """
}
