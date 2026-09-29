nextflow.enable.types = true

process CAT_FASTQ {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple(meta: Map, reads: List<Path>)

    stage:
    stageAs reads, 'input*/*'

    output:
    tuple(meta, files("*.merged.fastq.gz").toSorted { f -> f.name })

    topic:
    tuple(task.process, "cat", eval("cat --version 2>&1 | head -n 1 | sed 's/^.*coreutils) //; s/ .*\$//'")) >> 'versions'

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def compress = reads[0].name.endsWith('.gz') ? '' : '| gzip'
    if (meta.single_end) {
        if (reads.size() >= 1) {
            """
            cat ${reads.join(' ')} ${compress} > ${prefix}.merged.fastq.gz
            """
        } else {
            error("Could not find any FASTQ files to concatenate in the process input")
        }
    }
    else {
        if (reads.size() >= 2) {
            def read1 = reads.withIndex().findAll { _read, ix -> ix % 2 == 0 }.collect { read, _ix -> read }
            def read2 = reads.withIndex().findAll { _read, ix -> ix % 2 == 1 }.collect { read, _ix -> read }
            """
            cat ${read1.join(' ')} ${compress} > ${prefix}_1.merged.fastq.gz
            cat ${read2.join(' ')} ${compress} > ${prefix}_2.merged.fastq.gz
            """
        } else {
            error("Could not find any FASTQ file pairs to concatenate in the process input")
        }
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        if (reads.size() >= 1) {
            """
            echo '' | gzip > ${prefix}.merged.fastq.gz
            """
        } else {
            error("Could not find any FASTQ files to concatenate in the process input")
        }
    }
    else {
        if (reads.size() >= 2) {
            """
            echo '' | gzip > ${prefix}_1.merged.fastq.gz
            echo '' | gzip > ${prefix}_2.merged.fastq.gz
            """
        } else {
            error("Could not find any FASTQ file pairs to concatenate in the process input")
        }
    }
}
