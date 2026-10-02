nextflow.enable.types = true

include { GtfInput } from '../../types'

record EautilsGtf2bedResult {
    id:   String
    meta: Map
    bed:  Path
}

process EAUTILS_GTF2BED {
    tag "${sample.meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl:5.26.2' :
        'quay.io/biocontainers/perl:5.26.2' }"

    input:
    sample: GtfInput

    output:
    record(id: sample.id, meta: sample.meta, bed: file("${prefix}.bed")) as EautilsGtf2bedResult

    topic:
    file('versions.yml') >> 'versions'

    script:
    prefix = task.ext.prefix ?: "${sample.meta.id}"
    args   = task.ext.args ?: ''

    """
    echo $args
    """

    template 'gtf2bed.pl'

    stub:
    prefix = task.ext.prefix ?: "${sample.meta.id}"
    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl -e 'print substr(\$^V, 1)')
    END_VERSIONS
    """
}
