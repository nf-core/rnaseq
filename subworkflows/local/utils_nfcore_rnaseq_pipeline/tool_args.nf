//
// Tool arguments that follow from pipeline params. The workflows pass them to the processes that need
// them, so that a pipeline that includes this one controls them through its params record and a change
// only affects the tasks that use the value. Plain functions of values, so they can be shared by the
// processes of a tool.
//

// UMI-tools extract options
def umiExtractArgs(params) {
    return [
        params.umitools_extract_method ? "--extract-method=${params.umitools_extract_method}" : '',
        params.umitools_bc_pattern     ? "--bc-pattern='${params.umitools_bc_pattern}'" : '',
        params.umitools_bc_pattern2    ? "--bc-pattern2='${params.umitools_bc_pattern2}'" : '',
        params.umitools_umi_separator  ? "--umi-separator='${params.umitools_umi_separator}'" : ''
    ].join(' ').trim()
}
