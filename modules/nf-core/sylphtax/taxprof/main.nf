nextflow.enable.types = true

record SylphtaxTaxprofInput {
    id:          String
    meta:        Map
    profile_out: Path
}

process SYLPHTAX_TAXPROF {
    tag "${sample.meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/sylph-tax:1.9.0--pyhdfd78af_0'
        : 'quay.io/biocontainers/sylph-tax:1.9.0--pyhdfd78af_0'}"

    input:
    sample: SylphtaxTaxprofInput
    taxonomy: List<Path>

    output:
    record(id: sample.id, meta: sample.meta, taxprof_output: file('*.sylphmpa'))

    topic:
    tuple(task.process, 'sylph-tax', eval("sylph-tax --version 2>&1 | tail -1")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sample.meta.id}"

    """
    export SYLPH_TAXONOMY_CONFIG="/tmp/config.json"
    sylph-tax \\
        taxprof \\
        ${sample.profile_out} \\
        ${args} \\
        -t ${taxonomy.join(' ')}

    mv *.sylphmpa ${prefix}.sylphmpa
    """

    stub:
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    """
    touch ${prefix}.sylphmpa
    """
}
