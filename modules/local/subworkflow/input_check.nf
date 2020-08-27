/*
 * Check input samplesheet and get read channels
 */

include {
    SAMPLESHEET_CHECK;
    get_samplesheet_paths } from '../process/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet //   file: /path/to/samplesheet.csv
    seq_center  // string: sequencing center for read group
    options     //    map: options for samplesheet_check module

    main:
    SAMPLESHEET_CHECK ( samplesheet, options )
        .splitCsv ( header:true, sep:',' )
        .map { get_samplesheet_paths(it, seq_center) }
        .set { reads }

    emit:
    reads // channel: [ val(meta), [ reads ] ]
}
