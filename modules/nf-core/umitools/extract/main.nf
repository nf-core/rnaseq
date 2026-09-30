nextflow.enable.types = true

process UMITOOLS_EXTRACT {
    tag "$meta.id"
    label "process_single"
    label "process_long"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/32/32476f0107d72dbd2210a4e56b2873abde07300025cc11052680475509d2db81/data' :
        'community.wave.seqera.io/library/umi_tools_future_matplotlib_numpy_pruned:1ee668bafc8c9f81' }"

    input:
    tuple(meta: Map, reads: List<Path>)

    output:
    record(
        meta:  meta,
        reads: files('*.fastq.gz').toSorted { f -> f.name },
        log:   file('*.log')
    )

    topic:
    tuple(task.process, 'umitools', eval("umi_tools --version | sed -n '/version:/s/.*: //p'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        umi_tools \\
            extract \\
            -I ${[reads].flatten().join(' ')} \\
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
    def output_command = meta.single_end
        ? "echo '' | gzip > ${prefix}.umi_extract.fastq.gz"
        : "echo '' | gzip > ${prefix}.umi_extract_1.fastq.gz ;echo '' | gzip > ${prefix}.umi_extract_2.fastq.gz"
    """
    touch ${prefix}.umi_extract.log
    ${output_command}
    """
}
