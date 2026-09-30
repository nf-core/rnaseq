nextflow.enable.types = true

process HISAT2_ALIGN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/92/92f546b522b78597047a54609eefcf35a014de082df0e500bc2adf64d2733669/data'
        : 'community.wave.seqera.io/library/hisat2_samtools:a0c9b8ccf8116a89'}"

    input:
    tuple(meta: Map, reads: List<Path>)
    tuple(meta2: Map, index: Path)
    tuple(meta3: Map, splicesites: Path?)
    save_unaligned: Boolean

    output:
    record(
        id:       meta.id,
        meta:     meta,
        orig_bam: file('*.bam'),
        unmapped: files('*fastq.gz', optional: true).toSorted { f -> f.name },
        hisat2:   record(summary: file('*.log'))
    )

    topic:
    tuple(task.process, 'hisat2', eval("hisat2 --version | sed -n '1s/.*version //p'")) >> 'versions'
    tuple(task.process, 'samtools', eval("samtools --version | sed -n '1s/samtools //p'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def ss = splicesites ? "--known-splicesite-infile ${splicesites}" : ''
    def rg = args.contains("--rg-id") ? "" : "--rg-id ${prefix} --rg SM:${prefix}"
    if (meta.single_end) {
        def unaligned = save_unaligned ? "--un-gz ${prefix}.unmapped.fastq.gz" : ''
        """
        INDEX=`find -L ./ -name "*.1.ht2*" | sed 's/\\.1.ht2.*\$//'`
        hisat2 \\
            -x \$INDEX \\
            -U ${reads.join(' ')} \\
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
