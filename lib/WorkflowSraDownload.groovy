//
// This file holds functions specific to the workflow/sra_download.nf in the nf-core/rnaseq pipeline
//

class WorkflowSraDownload {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log, valid_params) {
        // Check minimal ENA fields are provided to download FastQ files
        def ena_metadata_fields = params.ena_metadata_fields ? params.ena_metadata_fields.split(',').collect{ it.trim().toLowerCase() } : valid_params['ena_metadata_fields']
        if (!ena_metadata_fields.containsAll(valid_params['ena_metadata_fields'])) {
            log.error "Invalid option: '${params.ena_metadata_fields}'. Minimally required fields for '--ena_metadata_fields': '${valid_params['ena_metadata_fields'].join(',')}'"
            System.exit(1)
        }
    }

    //
    // Print a warning after SRA download has completed
    //
    public static void sraDownloadWarn(log) {
        log.warn "=============================================================================\n" +
            "  Please double-check the samplesheet that has been auto-created using the\n" +
            "  public database ids provided via the '--public_data_ids' parameter.\n\n" +
            "  Public databases don't reliably hold information such as strandedness\n" +
            "  information.\n\n" +
            "  All of the sample metadata obtained from the ENA has been appended\n" +
            "  as additional columns to help you manually curate the samplesheet before\n" +
            "  you run the main workflow(s) in the pipeline.\n" +
            "==================================================================================="
    }
}
