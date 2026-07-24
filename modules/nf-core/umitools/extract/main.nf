process UMITOOLS_EXTRACT {
    tag "$meta.id"
    label "process_single"
    label "process_long"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/32/32476f0107d72dbd2210a4e56b2873abde07300025cc11052680475509d2db81/data' :
        'community.wave.seqera.io/library/umi_tools_future_matplotlib_numpy_pruned:1ee668bafc8c9f81' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    tuple val(meta), path("*.log")     , emit: log
    tuple val("${task.process}"), val('umitools'), eval("umi_tools --version | sed -n '/version:/s/.*: //p'"), emit: versions_umitools, topic: versions

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
    """
}
