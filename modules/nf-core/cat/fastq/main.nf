process CAT_FASTQ {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    meta    : Map
    reads   : List<Path>

    stage:
    stageAs "input*/*", reads

    output:
    file("*.merged.fastq.gz")

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def readList = reads.collect { it.toString() }
    if (meta.single_end) {
        if (readList.size >= 1) {
            """
            cat ${readList.join(' ')} > ${prefix}.merged.fastq.gz
            """
        }
    } else {
        if (readList.size >= 2) {
            def read1 = []
            def read2 = []
            readList.eachWithIndex{ v, ix -> ( ix & 1 ? read2 : read1 ) << v }
            """
            cat ${read1.join(' ')} > ${prefix}_1.merged.fastq.gz
            cat ${read2.join(' ')} > ${prefix}_2.merged.fastq.gz
            """
        }
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def readList = reads.collect { it.toString() }
    if (meta.single_end) {
        if (readList.size > 1) {
            """
            touch ${prefix}.merged.fastq.gz
            """
        }
    } else {
        if (readList.size > 2) {
            """
            touch ${prefix}_1.merged.fastq.gz
            touch ${prefix}_2.merged.fastq.gz
            """
        }
    }
}
