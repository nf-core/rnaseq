process BCFTOOLS_MPILEUP {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bcftools=1.16"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.16--hfe4b78e_1':
        'quay.io/biocontainers/bcftools:1.16--hfe4b78e_1' }"

    input:
    tuple val(meta), path(bam), path(intervals)
    path fasta
    val save_mpileup

    output:
    tuple val(meta), path("*vcf.gz")     , emit: vcf
    tuple val(meta), path("*vcf.gz.tbi") , emit: tbi
    tuple val(meta), path("*stats.txt")  , emit: stats
    tuple val(meta), path("*.mpileup.gz"), emit: mpileup, optional: true
    path  "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mpileup = save_mpileup ? "| tee ${prefix}.mpileup" : ""
    def bgzip_mpileup = save_mpileup ? "bgzip ${prefix}.mpileup" : ""
    def intervals = intervals ? "-T ${intervals}" : ""
    """
    echo "${meta.id}" > sample_name.list

    bcftools \\
        mpileup \\
        --fasta-ref $fasta \\
        $args \\
        $bam \\
        $intervals \\
        $mpileup \\
        | bcftools call --output-type v $args2 \\
        | bcftools reheader --samples sample_name.list \\
        | bcftools view --output-file ${prefix}.vcf.gz --output-type z $args3

    $bgzip_mpileup

    tabix -p vcf -f ${prefix}.vcf.gz

    bcftools stats ${prefix}.vcf.gz > ${prefix}.bcftools_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
