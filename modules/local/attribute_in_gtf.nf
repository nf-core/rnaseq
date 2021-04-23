/*
 * Check if attribute is available in column 9 of GTF file
 */
process ATTRIBUTE_IN_GTF {
    tag "$gtf"
    
    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    } else {
        container "biocontainers/biocontainers:v1.2.0_cv1"
    }

    input:
    path gtf
    val  attribute

    output:
    env  (ATTRIBUTE_EXISTS), emit: exists
    path ("*.txt")         , emit: txt
    
    script:
    """
    ATTRIBUTE_EXISTS=\$(awk -F "\t" '{ print \$9 }' $gtf | grep "$attribute " | wc -l | awk '{if(\$1 != 0) {print "true"} else print "false"}')
    echo \$ATTRIBUTE_EXISTS > attribute_exists.txt
    """
}