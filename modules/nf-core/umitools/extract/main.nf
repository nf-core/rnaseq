process UMITOOLS_EXTRACT {
    tag "$meta.id"
    label "process_single"
    label "process_long"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/umi_tools:1.1.5--py39hf95cd2a_0' :
        'biocontainers/umi_tools:1.1.5--py39hf95cd2a_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    tuple val(meta), path("*.log")     , emit: log
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        umi_tools \\
            extract \\
            -I $reads \\
            -S ${prefix}.umi_extract.fastq.gz \\
            $args \\
            > ${prefix}.umi_extract.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            umitools: \$( umi_tools --version | sed '/version:/!d; s/.*: //' )
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
        "${task.process}":
            umitools: \$( umi_tools --version | sed '/version:/!d; s/.*: //' )
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        output_command = "echo '' | gzip > ${prefix}.umi_extract.fastq.gz"
    } else {
        output_command = "echo '' | gzip > ${prefix}.umi_extract_1.fastq.gz ;"
        output_command += "echo '' | gzip > ${prefix}.umi_extract_2.fastq.gz"
    }
    """
    touch ${prefix}.umi_extract.log
    ${output_command}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        umitools: \$( umi_tools --version | sed '/version:/!d; s/.*: //' )
    END_VERSIONS
    """
}
