process UMITOOLS_EXTRACT {
    tag "$meta.id"
    label "process_low"

    conda (params.enable_conda ? "bioconda::umi_tools=1.1.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/umi_tools:1.1.2--py38h4a8c8d9_0"
    } else {
        container "quay.io/biocontainers/umi_tools:1.1.2--py38h4a8c8d9_0"
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    tuple val(meta), path("*.log")     , emit: log
    path  "versions.yml"               , emit: versions

    script:
    def prefix   = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def args     = task.ext.args ?: ''
    if (meta.single_end) {
        """
        umi_tools \\
            extract \\
            -I $reads \\
            -S ${prefix}.umi_extract.fastq.gz \\
            $args \\
            > ${prefix}.umi_extract.log

        cat <<-END_VERSIONS > versions.yml
        UMITOOLS_EXTRACT:
            umitools: \$(umi_tools --version 2>&1 | sed 's/^.*UMI-tools version://; s/ *\$//')
        END_VERSIONS
        """
    }  else {
        """
        umi_tools \\
            extract \\
            -I ${reads[0]} \\
            --read2-in=${reads[1]} \\
            -S ${prefix}.umi_extract_1.fastq.gz \\
            --read2-out=${prefix}.umi_extract_2.fastq.gz \\
            $args \\
            > ${prefix}.umi_extract.log

        cat <<-END_VERSIONS > versions.yml
        UMITOOLS_EXTRACT:
            umitools: \$(umi_tools --version 2>&1 | sed 's/^.*UMI-tools version://; s/ *\$//')
        END_VERSIONS
        """
    }
}
