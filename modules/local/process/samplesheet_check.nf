// Import generic module functions
include { saveFiles } from './functions'

/*
 * Reformat design file and check validity
 */
process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:'pipeline_info', publish_id:'') }

    conda (params.conda ? "${baseDir}/environment.yml" : null)

    input:
    path samplesheet
    val  options

    output:
    path '*.csv'

    script:  // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    """
    check_samplesheet.py $samplesheet samplesheet.valid.csv
    """
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def get_samplesheet_paths(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample
    meta.single_end   = row.single_end.toBoolean()
    meta.strandedness = row.strandedness

    def array = []
    if (meta.single_end) {
        array = [ meta, [ file(row.fastq_1, checkIfExists: true) ] ]
    } else {
        array = [ meta, [ file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true) ] ]
    }
    return array
}
