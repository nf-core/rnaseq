process TRIMGALORE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7e/7e44249e3fafe3d136ea726225551b51bca642387e16d9687b3e602207dedb20/data' :
        'community.wave.seqera.io/library/trim-galore:2.1.0--27e6376b8f6c1872'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*{3prime,5prime,trimmed,val}{,_1,_2}.fq.gz"), emit: reads
    tuple val(meta), path("*report.txt")                               , emit: log     , optional: true
    tuple val(meta), path("*unpaired{,_1,_2}.fq.gz")                   , emit: unpaired, optional: true
    tuple val(meta), path("*.html")                                    , emit: html    , optional: true
    tuple val(meta), path("*.zip")                                     , emit: zip     , optional: true
    tuple val("${task.process}"), val("trimgalore"), eval('trim_galore --version | grep -Eo "[0-9]+(\\.[0-9]+)+"'), topic: versions, emit: versions_trimgalore

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Calculate number of --cores for TrimGalore based on value of task.cpus
    // See: https://github.com/FelixKrueger/TrimGalore/blob/master/CHANGELOG.md#version-060-release-on-1-mar-2019
    // See: https://github.com/nf-core/atacseq/pull/65
    def cores = 1
    if (task.cpus) {
        cores = (task.cpus as int) - 4
        if (meta.single_end) {
            cores = (task.cpus as int) - 3
        }
        if (cores < 1) {
            cores = 1
        }
        if (cores > 8) {
            cores = 8
        }
    }

    // Added soft-links to original fastqs for consistent naming in MultiQC
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        def args_list = args.split("\\s(?=--)").toList()
        args_list.removeAll { arg -> arg.toLowerCase().contains('_r2 ') }
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -s ${reads} ${prefix}.fastq.gz
        trim_galore \\
            ${args_list.join(' ')} \\
            --cores ${cores} \\
            --gzip \\
            ${prefix}.fastq.gz
        """
    }
    else {
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
        trim_galore \\
            ${args} \\
            --cores ${cores} \\
            --paired \\
            --gzip \\
            ${prefix}_1.fastq.gz \\
            ${prefix}_2.fastq.gz
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        output_command = "echo '' | gzip > ${prefix}_trimmed.fq.gz ;"
        output_command += "touch ${prefix}.fastq.gz_trimming_report.txt"
    }
    else {
        output_command = "echo '' | gzip > ${prefix}_1_trimmed.fq.gz ;"
        output_command += "touch ${prefix}_1.fastq.gz_trimming_report.txt ;"
        output_command += "echo '' | gzip > ${prefix}_2_trimmed.fq.gz ;"
        output_command += "touch ${prefix}_2.fastq.gz_trimming_report.txt"
    }
    """
    ${output_command}
    """
}
