nextflow.enable.types = true

record CustomMultiqccustombiotypeInput {
    id:     String
    meta:   Map
    counts: Path
}

record CustomMultiqccustombiotypeResult {
    id:   String
    meta: Map
    tsv:  Path
    rrna: Path
}

process CUSTOM_MULTIQCCUSTOMBIOTYPE {
    tag "${sample.meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12.12' :
        'quay.io/biocontainers/python:3.12.12' }"

    input:
    sample: CustomMultiqccustombiotypeInput
    header: Path

    output:
    record(
        id:   sample.id,
        meta: sample.meta,
        tsv:  file('*biotype_counts_mqc.tsv'),
        rrna: file('*biotype_counts_rrna_mqc.tsv')
    ) as CustomMultiqccustombiotypeResult

    topic:
    file('versions.yml') >> 'versions'

    script:
    template 'mqc_features_stat.py'

    stub:
    def prefix = task.ext.prefix ?: "${sample.meta.id}"
    """
    touch ${prefix}.biotype_counts_mqc.tsv
    touch ${prefix}.biotype_counts_rrna_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
