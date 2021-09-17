//
// Convert tab file from STAR to junc format
//

// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process TAB_TO_JUNC {
    input:
    path tab

    output:
    tuple val(meta), path('*.junc')       , emit: junc

    """
    cp $tab > ${tab}.junc
    """
}
