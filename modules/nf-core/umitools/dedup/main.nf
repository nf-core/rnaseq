process UMITOOLS_DEDUP {
    tag "$meta.id"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/umi_tools:1.1.5--py39hf95cd2a_0' :
        'biocontainers/umi_tools:1.1.5--py39hf95cd2a_0' }"

    input:
    meta                : Map
    bam                 : Path
    bai                 : Path
    get_output_stats    : boolean

    output:
    bam                     : Path = file("${prefix}.bam")
    log                     : Path = file("*.log")
    tsv_edit_distance       : Path? = file("*_edit_distance.tsv")
    tsv_per_umi             : Path? = file("*_per_umi.tsv")
    tsv_umi_per_position    : Path? = file("*_per_position.tsv")

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def paired = meta.single_end ? "" : "--paired"
    stats = get_output_stats ? "--output-stats ${prefix}" : ""
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    if (!(args ==~ /.*--random-seed.*/)) {args += " --random-seed=100"}
    """
    PYTHONHASHSEED=0 umi_tools \\
        dedup \\
        -I $bam \\
        -S ${prefix}.bam \\
        -L ${prefix}.log \\
        $stats \\
        $paired \\
        $args
    """

    stub:
    """
    touch ${prefix}.bam
    touch ${prefix}.log
    touch ${prefix}_edit_distance.tsv
    touch ${prefix}_per_umi.tsv
    touch ${prefix}_per_position.tsv
    """
}
