process UMITOOLS_DEDUP {
    tag "$meta.id"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/32/32476f0107d72dbd2210a4e56b2873abde07300025cc11052680475509d2db81/data' :
        'community.wave.seqera.io/library/umi_tools_future_matplotlib_numpy_pruned:1ee668bafc8c9f81' }"

    input:
    tuple val(meta), path(bam), path(bai)
    val get_output_stats

    output:
    tuple val(meta), path("${prefix}.bam")     , emit: bam
    tuple val(meta), path("*.log")             , emit: log
    tuple val(meta), path("*edit_distance.tsv"), optional:true, emit: tsv_edit_distance
    tuple val(meta), path("*per_umi.tsv")      , optional:true, emit: tsv_per_umi
    tuple val(meta), path("*per_position.tsv") , optional:true, emit: tsv_umi_per_position
    tuple val("${task.process}"), val('umitools'), eval("umi_tools --version | sed 's/UMI-tools version: //'"), emit: versions_umitools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def paired = meta.single_end ? "" : "--paired"
    stats = get_output_stats ? "--output-stats ${prefix}" : ""
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    if (!(args ==~ /.*--random-seed.*/)) {args += " --random-seed=100"}
    """
    #Prevent matplotlib from using /tmp
    mkdir .tmp && chmod 777 .tmp

    MPLCONFIGDIR=.tmp TMPDIR=.tmp PYTHONHASHSEED=0 umi_tools \\
        dedup \\
        -I $bam \\
        -S ${prefix}.bam \\
        -L ${prefix}.log \\
        $stats \\
        $paired \\
        $args
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.log
    touch ${prefix}_edit_distance.tsv
    touch ${prefix}_per_umi.tsv
    touch ${prefix}_per_position.tsv
    """
}
