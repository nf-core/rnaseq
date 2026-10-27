process HISAT2_ALIGN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/92/92f546b522b78597047a54609eefcf35a014de082df0e500bc2adf64d2733669/data'
        : 'community.wave.seqera.io/library/hisat2_samtools:a0c9b8ccf8116a89'}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    tuple val(meta3), path(splicesites)
    val save_unaligned

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.log"), emit: summary
    tuple val(meta), path("*fastq.gz"), optional: true, emit: fastq
    tuple val("${task.process}"), val('hisat2'), eval("hisat2 --version | sed -n '1s/.*version //p'"), emit: versions_hisat2, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | sed -n '1s/samtools //p'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def ss = "${splicesites}" ? "--known-splicesite-infile ${splicesites}" : ''
    def rg = args.contains("--rg-id") ? "" : "--rg-id ${prefix} --rg SM:${prefix}"
    if (meta.single_end) {
        def unaligned = save_unaligned ? "--un-gz ${prefix}.unmapped.fastq.gz" : ''
        """
        INDEX=`find -L ./ -name "*.1.ht2*" | sed 's/\\.1.ht2.*\$//'`
        hisat2 \\
            -x \$INDEX \\
            -U ${reads} \\
            ${ss} \\
            --summary-file ${prefix}.hisat2.summary.log \\
            --threads ${task.cpus} \\
            ${rg} \\
            ${unaligned} \\
            ${args} \\
            | samtools view -bS -F 4 -F 256 - > ${prefix}.bam
        """
    }
    else {
        def unaligned = save_unaligned ? "--un-conc-gz ${prefix}.unmapped.fastq.gz" : ''
        """
        INDEX=`find -L ./ -name "*.1.ht2*" | sed 's/\\.1.ht2.*\$//'`
        hisat2 \\
            -x \$INDEX \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            ${ss} \\
            --summary-file ${prefix}.hisat2.summary.log \\
            --threads ${task.cpus} \\
            ${rg} \\
            ${unaligned} \\
            --no-mixed \\
            --no-discordant \\
            ${args} \\
            | samtools view -bS -F 4 -F 8 -F 256 - > ${prefix}.bam

        if [ -f ${prefix}.unmapped.fastq.1.gz ]; then
            mv ${prefix}.unmapped.fastq.1.gz ${prefix}.unmapped_1.fastq.gz
        fi
        if [ -f ${prefix}.unmapped.fastq.2.gz ]; then
            mv ${prefix}.unmapped.fastq.2.gz ${prefix}.unmapped_2.fastq.gz
        fi
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def unaligned = save_unaligned ? "echo '' | gzip >  ${prefix}.unmapped_1.fastq.gz \n echo '' | gzip >  ${prefix}.unmapped_2.fastq.gz" : ''
    """
    ${unaligned}

    touch ${prefix}.hisat2.summary.log
    touch ${prefix}.bam
    """
}
