// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

/*
 * Convert GTF file to BED format
 */
process GTF2BED {
    tag "$gtf"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'genome', publish_id:'') }

    container "quay.io/biocontainers/perl:5.26.2"
    //container "https://depot.galaxyproject.org/singularity/perl:5.26.2"

    conda (params.conda ? "conda-forge::perl=5.26.2" : null)

    input:
    path gtf
    
    output:
    path '*.bed'

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    """
    gtf2bed $gtf > ${gtf.baseName}.bed
    """
}
