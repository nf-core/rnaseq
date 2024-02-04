process RSEQC_INNERDISTANCE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::rseqc=5.0.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rseqc:5.0.3--py39hf95cd2a_0' :
        'biocontainers/rseqc:5.0.3--py39hf95cd2a_0' }"

    input:
    tuple val(meta), path(bam)
    path  bed

    output:
    tuple val(meta), path("*distance.txt"), optional:true, emit: distance
    tuple val(meta), path("*freq.txt")    , optional:true, emit: freq
    tuple val(meta), path("*mean.txt")    , optional:true, emit: mean
    tuple val(meta), path("*.pdf")        , optional:true, emit: pdf
    tuple val(meta), path("*.r")          , optional:true, emit: rscript
    path  "versions.yml"                  , emit: versions

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

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            rseqc: \$(inner_distance.py --version | sed -e "s/inner_distance.py //g")
        END_VERSIONS
        """
    } else {
        """
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            rseqc: \$(inner_distance.py --version | sed -e "s/inner_distance.py //g")
        END_VERSIONS
        """
    }
}
