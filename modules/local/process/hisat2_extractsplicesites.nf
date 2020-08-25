// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process HISAT2_EXTRACTSPLICESITES {
    tag "$gtf"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container "quay.io/biocontainers/hisat2:2.2.0--py37hfa133b6_4"
    //container "https://depot.galaxyproject.org/singularity/hisat2:2.2.0--py37hfa133b6_4"
    
    conda (params.conda ? "bioconda::hisat2=2.2.0" : null)

    input:
    path gtf
    val options

    output:
    path "*.splice_sites.txt", emit: txt
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    hisat2_extract_splice_sites.py $gtf > ${gtf.baseName}.splice_sites.txt
    echo \$(hisat2 --version 2>&1) | sed 's/^.*version //; s/64.*\$//' > ${software}.version.txt
    """
}
