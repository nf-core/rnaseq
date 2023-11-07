process GTF_FOR_STRINGTIE {
    tag "$gtf"
    label 'process_low'

    conda "conda-forge::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'biocontainers/gawk:5.1.0' }"

    input:
    path gtf

    output:
    path '*_for_stringtie.gtf' , emit: gtf
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    """
    # Remove gene entries from a GTF, if they don't contain a transcript_id
    # attribute, which is common in Ensembl GTFs
    awk '\$3 == "gene" { print; next } /transcript_id/ && \$0 ~ /transcript_id "[^"]+"/' $gtf > ${gtf}_for_stringtie.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(echo \$(bash --version | grep -Eo 'version [[:alnum:].]+' | sed 's/version //'))
    END_VERSIONS
    """
}
