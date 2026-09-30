nextflow.enable.types = true

process EAUTILS_GTF2BED {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl:5.26.2' :
        'quay.io/biocontainers/perl:5.26.2' }"

    input:
    record(id: String, meta: Map, gtf: Path)

    output:
    record(id: id, meta: meta, bed: file("${prefix}.bed"))

    topic:
    file('versions.yml') >> 'versions'

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    args   = task.ext.args ?: ''

    """
    echo $args
    """

    template 'gtf2bed.pl'

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl -e 'print substr(\$^V, 1)')
    END_VERSIONS
    """
}
