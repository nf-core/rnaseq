process RSEQC_INNERDISTANCE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6f/6f44b7933e2c2b1a340dc9485869974eb032d34e81af83716eb381964ee3e5e7/data' :
        'community.wave.seqera.io/library/rseqc_r-base:2e29d2dfda9cef15' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  bed

    output:
    tuple val(meta), path("*distance.txt"), optional:true, emit: distance
    tuple val(meta), path("*freq.txt")    , optional:true, emit: freq
    tuple val(meta), path("*mean.txt")    , optional:true, emit: mean
    tuple val(meta), path("*.pdf")        , optional:true, emit: pdf
    tuple val(meta), path("*.r")          , optional:true, emit: rscript
    tuple val("${task.process}"), val('rseqc'), eval('inner_distance.py --version | sed "s/inner_distance.py //"'), emit: versions_rseqc, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (!meta.single_end) {
        """
        inner_distance.py \\
            -i $bam \\
            -r $bed \\
            -o $prefix \\
            $args \\
            > stdout.txt
        head -n 2 stdout.txt > ${prefix}.inner_distance_mean.txt
        """
    } else {
        """
        echo "inner_distance.py doesn't support single-end data" > ${prefix}.inner_distance.txt
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.inner_distance.txt
    touch ${prefix}.inner_distance_freq.txt
    touch ${prefix}.inner_distance_mean.txt
    touch ${prefix}.inner_distance_plot.pdf
    touch ${prefix}.inner_distance_plot.r
    """
}
