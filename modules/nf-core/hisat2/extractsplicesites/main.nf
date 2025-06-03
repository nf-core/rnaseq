process HISAT2_EXTRACTSPLICESITES {
    tag "$gtf"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hisat2:2.2.1--h1b792b2_3' :
        'biocontainers/hisat2:2.2.1--h1b792b2_3' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.splice_sites.txt"), emit: txt
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    hisat2_extract_splice_sites.py $gtf > ${gtf.baseName}.splice_sites.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2: \$(hisat2 --version | grep -o 'version [^ ]*' | cut -d ' ' -f 2)
    END_VERSIONS
    """

    stub:
    """
    touch ${gtf.baseName}.splice_sites.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2: \$(hisat2 --version | grep -o 'version [^ ]*' | cut -d ' ' -f 2)
    END_VERSIONS
    """
}
