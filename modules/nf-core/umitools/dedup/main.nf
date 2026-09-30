nextflow.enable.types = true

include { BamBaiInput } from '../../types'

process UMITOOLS_DEDUP {
    tag "${sample.meta.id}"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/32/32476f0107d72dbd2210a4e56b2873abde07300025cc11052680475509d2db81/data' :
        'community.wave.seqera.io/library/umi_tools_future_matplotlib_numpy_pruned:1ee668bafc8c9f81' }"

    input:
    sample: BamBaiInput
    get_output_stats: Boolean

    output:
    record(
        id:                   sample.id,
        meta:                 sample.meta,
        bam:                  file("${prefix}.bam"),
        log:                  file('*.log'),
        tsv_edit_distance:    file('*edit_distance.tsv', optional: true),
        tsv_per_umi:          file('*per_umi.tsv', optional: true),
        tsv_umi_per_position: file('*per_position.tsv', optional: true)
    )

    topic:
    tuple(task.process, 'umitools', eval("umi_tools --version | sed 's/UMI-tools version: //'")) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${sample.meta.id}"
    def paired = sample.meta.single_end ? "" : "--paired"
    stats = get_output_stats ? "--output-stats ${prefix}" : ""
    if ("${sample.bam}" == "${prefix}.bam") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    if (!(args ==~ /.*--random-seed.*/)) {args += " --random-seed=100"}
    """
    #Prevent matplotlib from using /tmp
    mkdir .tmp && chmod 777 .tmp

    MPLCONFIGDIR=.tmp TMPDIR=.tmp PYTHONHASHSEED=0 umi_tools \\
        dedup \\
        -I ${sample.bam} \\
        -S ${prefix}.bam \\
        -L ${prefix}.log \\
        $stats \\
        $paired \\
        $args
    """

    stub:
    prefix = task.ext.prefix ?: "${sample.meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.log
    touch ${prefix}_edit_distance.tsv
    touch ${prefix}_per_umi.tsv
    touch ${prefix}_per_position.tsv
    """
}
