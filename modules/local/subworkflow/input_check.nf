/*
 * Check input samplesheet and get read channels
 */

include { SAMPLESHEET_CHECK;
          get_samplesheet_paths } from '../process/samplesheet_check'

workflow INPUT_CHECK {
    take:
    ch_input                  //   file: /path/to/samplesheet.csv
    seq_center                // string: sequencing center for read group
    samplesheet_check_options //    map: options for check_samplesheet module

    main:
    SAMPLESHEET_CHECK (ch_input, samplesheet_check_options)
        .splitCsv(header:true, sep:',')
        .map { get_samplesheet_paths(it, seq_center) }
        .set { ch_reads }

    emit:
    reads = ch_reads // channel: [ val(meta), [ reads ] ]
}
