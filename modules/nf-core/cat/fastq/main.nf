process CAT_FASTQ {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple val(meta), path(reads, stageAs: "input*/*")

    output:
    tuple val(meta), path("*.merged.fastq.gz"), emit: reads
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def readList = reads instanceof List ? reads.collect { it.toString() } : [reads.toString()]
    if (meta.single_end) {
        if (readList.size >= 1) {
            """
            cat ${readList.join(' ')} > ${prefix}.merged.fastq.gz

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
            END_VERSIONS
            """
        } else {
            error("Could not find any FASTQ files to concatenate in the process input")
        }
    }
    else {
        if (readList.size >= 2) {
            def read1 = []
            def read2 = []
            readList.eachWithIndex { v, ix -> (ix & 1 ? read2 : read1) << v }
            """
            cat ${read1.join(' ')} > ${prefix}_1.merged.fastq.gz
            cat ${read2.join(' ')} > ${prefix}_2.merged.fastq.gz

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
            END_VERSIONS
            """
        } else {
            error("Could not find any FASTQ file pairs to concatenate in the process input")
        }
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def readList = reads instanceof List ? reads.collect { it.toString() } : [reads.toString()]
    if (meta.single_end) {
        if (readList.size >= 1) {
            """
            echo '' | gzip > ${prefix}.merged.fastq.gz

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
            END_VERSIONS
            """
        } else {
            error("Could not find any FASTQ files to concatenate in the process input")
        }
    }
    else {
        if (readList.size >= 2) {
            """
            echo '' | gzip > ${prefix}_1.merged.fastq.gz
            echo '' | gzip > ${prefix}_2.merged.fastq.gz

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
            END_VERSIONS
            """
        } else {
            error("Could not find any FASTQ file pairs to concatenate in the process input")
        }
    }
}
